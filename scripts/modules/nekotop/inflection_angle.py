import numpy as np


# Function definition for the angle computation
def compute_angle(center, point, axis=np.array([1.0, 0.0, 0.0])):
    """
    Compute the angle between the center and the point
    """

    # Compute the angle
    angle = np.arctan2(point[1] - center[1], point[0] - center[0])

    return angle

def compute_radial_flow(angle, field) -> float:
    """
    Compute the radial flow of the field
    """

    # Compute the flow
    flow = -field[0] * np.sin(angle) + field[1] * np.cos(angle)

    return flow

def compute_inflection_angle(i0, direction, flow, angles) -> float:

    i_p = i0
    i_n = i_p + direction

    while flow[i_p] < 0 or flow[i_n] > 0:
        if flow[i_n] < 0:
            s = - direction
        else:
            s = + direction

        i_p = (i_p + s) % len(angles)
        i_n = (i_n + s) % len(angles)

    # Compute the angle of the inflection point
    angle_p = angles[i_p]
    angle_n = angles[i_n]
    flow_p = flow[i_p]
    flow_n = flow[i_n]

    return angle_n + (angle_p - angle_n) * flow_n / (flow_n - flow_p)

def track_inflection_point(center, points, fields, times):
    """
    Track the inflection point of the flow
    """

    N_points = points.shape[0]
    N_times = times.shape[0]

    angle_inflection = np.zeros([len(times), 3])
    angles = np.zeros(points.shape[0])
    for i_boundary in range(points.shape[0]):
        point = points[i_boundary, :]
        angles[i_boundary] = compute_angle(center, point)

    idx = np.argsort(angles)

    for ti in range(N_times):
        t = times[ti]

        flow = np.zeros(points.shape[0])
        for i_boundary in range(points.shape[0]):
            angle = angles[i_boundary]
            field = fields[i_boundary + ti * N_points, :]
            flow[i_boundary] = compute_radial_flow(angle, field)

        i0 = int(3 / 4 * len(idx))
        a0 = compute_inflection_angle(i0, 1, flow[idx], angles[idx])

        i0 = int(2 / 4 * len(idx))
        a1 = compute_inflection_angle(i0, -1, flow[idx], angles[idx])

        i0 = int(1 / 4 * len(idx))
        a2 = compute_inflection_angle(i0, 1, flow[idx], angles[idx])

        # The central inflection point cannot be outside the other two.
        if a0 < a1 or a1 < a2:
            a1 = 0.5 * (a0 + a2)

        angle_inflection[ti, :] = [a0, a1, a2]

    return angle_inflection

