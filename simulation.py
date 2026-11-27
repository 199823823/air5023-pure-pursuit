"""Fixed-lookahead Pure Pursuit simulation for the AIR6023 group project."""

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Settings:
    velocity_mps: float = 1.0
    wheelbase_m: float = 1.0
    time_step_s: float = 0.05
    steering_limit_deg: float = 45.0
    initial_lateral_offset_m: float = 0.5
    completion_tolerance_m: float = 0.25
    max_time_s: float = 60.0


SETTINGS = Settings()
PATH_NAMES = {"straight": "Straight", "circle": "Circular", "s_curve": "S-shaped"}
LOOKAHEAD_VALUES = (0.5, 1.0, 2.0, 3.0, 4.0)
GEOMETRY_TOL = 1e-9


class ReferencePath:
    """An ordered polyline with cumulative chord lengths."""

    def __init__(self, points, start_heading, closed=False):
        self.points = points
        self.start_heading = start_heading
        self.closed = closed
        self.starts = points[:-1]
        self.vectors = np.diff(points, axis=0)
        self.length_sq = np.sum(self.vectors ** 2, axis=1)
        self.lengths = np.sqrt(self.length_sq)
        self.arc = np.r_[0.0, np.cumsum(self.lengths)]
        self.length = self.arc[-1]

    def project(self, position):
        """Return the minimum distance, closest segment index and arc position."""
        offsets = position - self.starts
        fractions = np.clip(
            np.sum(offsets * self.vectors, axis=1) / self.length_sq, 0.0, 1.0
        )
        residuals = offsets - fractions[:, None] * self.vectors
        distances_sq = np.sum(residuals ** 2, axis=1)
        index = int(np.argmin(distances_sq))
        distance = np.sqrt(distances_sq[index])
        arc_position = self.arc[index] + fractions[index] * self.lengths[index]
        return distance, index, arc_position

    def point_at(self, arc_position):
        if self.closed:
            arc_position %= self.length
        return np.array(
            [np.interp(arc_position, self.arc, self.points[:, axis]) for axis in range(2)]
        )

    def target(self, position, segment, arc_position, radius):
        """Find the first forward intersection on the exiting side of the disk."""
        count = len(self.starts)
        segments = np.arange(segment, count)
        lap_offsets = np.zeros(len(segments))
        if self.closed:
            segments = np.r_[segments, np.arange(segment + 1)]
            lap_offsets = np.r_[lap_offsets, np.full(segment + 1, self.length)]

        vectors = self.vectors[segments]
        offsets = self.starts[segments] - position
        length_sq = self.length_sq[segments]
        # Intersect each segment's supporting line with the lookahead circle.
        center_fraction = -np.sum(offsets * vectors, axis=1) / length_sq
        normals = offsets + center_fraction[:, None] * vectors
        room = radius ** 2 - np.sum(normals ** 2, axis=1)
        exits = center_fraction + np.sqrt(np.maximum(room, 0.0) / length_sq)
        target_arc = self.arc[segments] + lap_offsets + exits * self.lengths[segments]
        valid = (
            (room >= -GEOMETRY_TOL)
            & (exits >= -GEOMETRY_TOL)
            & (exits <= 1.0 + GEOMETRY_TOL)
            & (target_arc >= arc_position - GEOMETRY_TOL)
            & (target_arc <= arc_position + self.length + GEOMETRY_TOL)
        )
        matches = np.flatnonzero(valid)
        if matches.size:
            first = matches[0]
            return self.starts[segments[first]] + np.clip(exits[first], 0, 1) * vectors[first]

        if not self.closed and np.linalg.norm(self.points[-1] - position) <= radius:
            return self.points[-1]
        # Advance along the path during recovery if it misses the disk.
        return self.point_at(arc_position + radius)


def make_path(name):
    if name == "straight":
        x = np.linspace(0.0, 32.0, 801)
        return ReferencePath(np.column_stack((x, np.zeros_like(x))), 0.0)
    if name == "circle":
        phi = np.linspace(0.0, 2.0 * np.pi, 801)
        points = np.column_stack((5.0 * np.sin(phi), 5.0 * (1.0 - np.cos(phi))))
        points[-1] = points[0]
        return ReferencePath(points, 0.0, closed=True)
    if name == "s_curve":
        x = np.linspace(0.0, 32.0, 801)
        frequency = 2.0 * np.pi / 16.0
        return ReferencePath(
            np.column_stack((x, 2.0 * np.sin(frequency * x))),
            np.arctan(2.0 * frequency),
        )
    raise ValueError(f"Unknown reference path: {name}")


@dataclass
class SimulationResult:
    lookahead_m: float
    time_s: np.ndarray
    poses: np.ndarray
    tracking_error_m: np.ndarray
    progress_m: np.ndarray
    targets: np.ndarray
    target_distance_m: np.ndarray
    steering_deg: np.ndarray
    completed: bool

    def metrics(self, path_name):
        # Commands occupy integration intervals. The last state has no command.
        steering_changes = np.abs(np.diff(self.steering_deg))
        variation_rate = (
            float(steering_changes.sum() / (self.time_s[len(self.steering_deg) - 1]))
            if steering_changes.size else 0.0
        )
        return {
            "path": path_name,
            "lookahead_m": self.lookahead_m,
            "rmse_m": float(np.sqrt(np.mean(self.tracking_error_m ** 2))),
            "max_tracking_error_m": float(self.tracking_error_m.max()),
            "steering_variation_rate_deg_per_s": variation_rate,
            "path_completed": self.completed,
            "elapsed_time_s": float(self.time_s[-1]),
        }


def simulate(path, lookahead_m, settings=SETTINGS):
    theta = path.start_heading
    left_normal = np.array([-np.sin(theta), np.cos(theta)])
    x, y = path.points[0] + settings.initial_lateral_offset_m * left_normal
    dt = settings.time_step_s
    step_count = int(np.floor(settings.max_time_s / dt))
    steering_limit = np.deg2rad(settings.steering_limit_deg)
    travel_per_step = settings.velocity_mps * dt
    previous_arc = 0.0
    completed = False
    poses, errors, progress_log, targets, distances, steering_log = [], [], [], [], [], []

    for step in range(step_count + 1):
        position = np.array([x, y])
        error, segment, arc_position = path.project(position)
        progress = arc_position
        tolerance = settings.completion_tolerance_m
        if path.closed:
            # Choose the lap copy nearest the previous projection across the seam.
            laps = np.rint((previous_arc - arc_position) / path.length)
            progress += laps * path.length
            completed = progress >= path.length and error <= tolerance
        else:
            endpoint_distance = np.linalg.norm(position - path.points[-1])
            completed = path.length - progress <= tolerance and endpoint_distance <= tolerance
        previous_arc = progress

        poses.append((x, y, theta))
        errors.append(error)
        progress_log.append(progress)
        if completed or step == step_count:
            break

        target = path.target(position, segment, arc_position, lookahead_m)
        dx, dy = target - position
        target_distance = np.hypot(dx, dy)
        alpha = np.arctan2(dy, dx) - theta
        curvature = 2.0 * np.sin(alpha) / target_distance if target_distance > GEOMETRY_TOL else 0.0
        steering = np.clip(np.arctan(settings.wheelbase_m * curvature), -steering_limit, steering_limit)
        targets.append(target)
        distances.append(target_distance)
        steering_log.append(steering)

        x += travel_per_step * np.cos(theta)
        y += travel_per_step * np.sin(theta)
        theta += travel_per_step * np.tan(steering) / settings.wheelbase_m
        theta = (theta + np.pi) % (2.0 * np.pi) - np.pi

    return SimulationResult(
        lookahead_m, np.arange(len(poses)) * dt, np.asarray(poses),
        np.asarray(errors), np.asarray(progress_log), np.asarray(targets).reshape(-1, 2),
        np.asarray(distances), np.rad2deg(steering_log), completed,
    )


def write_trajectory(file, result):
    with file.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "time_s", "x_m", "y_m", "heading_rad", "tracking_error_m", "progress_m",
            "target_x_m", "target_y_m", "target_distance_m", "steering_deg",
        ])
        for i, time in enumerate(result.time_s):
            state = [time, *result.poses[i], result.tracking_error_m[i], result.progress_m[i]]
            control = ["", "", "", ""]
            if i < len(result.steering_deg):
                control = [*result.targets[i], result.target_distance_m[i], result.steering_deg[i]]
            writer.writerow([*state, *control])


def plot_results(file, path_name, path, results):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig = plt.figure(figsize=(11, 6), layout="constrained")
    grid = fig.add_gridspec(2, 2)
    if path.closed:
        trajectory = fig.add_subplot(grid[:, 0])
        error = fig.add_subplot(grid[0, 1])
        steering = fig.add_subplot(grid[1, 1], sharex=error)
        error.tick_params(labelbottom=False)
    else:
        trajectory = fig.add_subplot(grid[0, :])
        error = fig.add_subplot(grid[1, 0])
        steering = fig.add_subplot(grid[1, 1], sharex=error)
        error.set_xlabel("Time (s)")
    trajectory.plot(*path.points.T, color="0.25", linestyle="--", linewidth=1.5, label="Reference")
    for result in results:
        label = f"{result.lookahead_m:.1f} m"
        line, = trajectory.plot(*result.poses[:, :2].T, label=label)
        error.plot(result.time_s, result.tracking_error_m, color=line.get_color())
        if result.steering_deg.size:
            steering.step(
                result.time_s, np.r_[result.steering_deg, result.steering_deg[-1]],
                where="post", color=line.get_color(),
            )
    trajectory.set(xlabel="x (m)", ylabel="y (m)")
    lower, upper = trajectory.get_ylim()
    if upper - lower < 2.0:
        center = (lower + upper) / 2.0
        trajectory.set_ylim(center - 1.0, center + 1.0)
    trajectory.set_aspect("equal", adjustable="box")
    trajectory.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01),
        ncols=min(3 if path.closed else 6, len(results) + 1), frameon=False,
    )
    error.set(title="Tracking error", ylabel="Distance (m)", ylim=(0.0, None))
    steering.set(title="Applied steering", xlabel="Time (s)", ylabel="Angle (deg)")
    for axis in (trajectory, error, steering):
        axis.grid(alpha=0.25)
    fig.suptitle(f"Pure Pursuit | {PATH_NAMES[path_name]} path")
    fig.savefig(file, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", choices=PATH_NAMES, default="s_curve")
    parser.add_argument(
        "--lookahead", nargs="+", type=float, choices=LOOKAHEAD_VALUES, default=[2.0],
        help="One or more lookahead radii in metres (default: 2.0).",
    )
    parser.add_argument(
        "--output", type=Path, help="Output directory (default: outputs/<path>).",
    )
    args = parser.parse_args()
    output = args.output or Path("outputs") / args.path
    path = make_path(args.path)
    results = [simulate(path, value) for value in dict.fromkeys(args.lookahead)]

    output.mkdir(parents=True, exist_ok=True)
    summaries = [result.metrics(args.path) for result in results]
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    for result in results:
        write_trajectory(output / f"trajectory_{result.lookahead_m:.1f}m.csv", result)
    plot_results(output / "tracking.png", args.path, path, results)
    for summary in summaries:
        print(
            f"{summary['path']}, {summary['lookahead_m']:.1f} m: "
            f"completed={summary['path_completed']}, RMSE={summary['rmse_m']:.4f} m, "
            f"time={summary['elapsed_time_s']:.2f} s"
        )
    print(f"Results written to {output.resolve()}")


if __name__ == "__main__":
    main()
