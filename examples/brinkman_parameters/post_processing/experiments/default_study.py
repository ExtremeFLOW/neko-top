import numpy as np
import matplotlib.pyplot as plt
def plot_study(experiment_tabulated, x_axis_variable, line_axis_variable, out_filename):
# its late at night and im tired, so here are the todo's:
# - allow for an optional input to hold an aditional variable constant, ie, Re
# - include the separation angle observable
# - handle the axis names etc
# - make the actual plotting neat and tidy!
    # for now I'm hardcoding a list of observables
    lift_obs = {
        "type"      : "force",
        "label"     : "Lift",
        "circle"    : "050",
        "component" : "fy_tot",
        "linestyle" : "-",
        "axis"      : 1,
        }
    drag_obs = {
        "type"      : "force",
        "label"     : "Drag",
        "circle"    : "050",
        "component" : "fx_tot",
        "linestyle" : "--",
        "axis"      : 2,
        }
    observable_list = [lift_obs, drag_obs]
    # get the plot read
    fig, ax = plt.subplots()
    # to be more general we should check how many axis are in the observable lsit
    # I'm no good with python, I want something along the lines of plotyy in matlab
    ax1 = ax.twinx()
    ax2 = ax.twinx()
    # keep the colours
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    # ok, we're going to look for the indipendent variables
    method = []
    mesh = []
    Re = []
    chi = []
    radius = []
    rmax = []
    for case in experiment_tabulated:
        method.append(case["method"])
        mesh.append(case["mesh"])
        Re.append(case["Re"])
        chi.append(case["chi"])
        radius.append(case["radius"])
        rmax.append(case["rmax"])
    unique_method = list(set(method))
    if None in unique_method:
        unique_method.remove(None)
    unique_mesh = list(set(mesh))
    if None in unique_mesh:
        unique_mesh.remove(None)
    unique_Re = list(set(Re))
    if None in unique_Re:
        unique_Re.remove(None)
    unique_chi = list(set(chi))
    if None in unique_chi:
        unique_chi.remove(None)
    unique_radius = list(set(radius))
    if None in unique_radius:
        unique_radius.remove(None)
    unique_rmax = list(set(rmax))
    if None in unique_rmax:
        unique_rmax.remove(None)
    # This is a bit tricky, but in my opinion we can compare at most 2 indipendent
    # variables in a 2D plot, otherwise we need to move to 3D.
    # our options are:
    #  - use the x axis
    #  - use multiple lines

    # Tim, I know you wanted to the put the radius and rmax both on the x axis
    # But that's a bit strange to me...
    # If you still want to do it I guess we can, but it seems a bit inconsistent
    variable_list = ["method", "mesh", "Re", "chi", "radius", "rmax"]           
    unique_list = [unique_method, unique_mesh, unique_Re, unique_chi,           
        unique_radius, unique_rmax]  
    ind = variable_list.index(x_axis_variable)
    x_axis_list = unique_list[ind]
    ind = variable_list.index(line_axis_variable)
    line_list = unique_list[ind]
    for j, line in enumerate(line_list):
        observables = np.empty((len(x_axis_list), len(line_list)))
        for i, x_pt in enumerate(x_axis_list):
            for c, case in enumerate(experiment_tabulated):
                if line == case[line_axis_variable] and x_pt == case[x_axis_variable]:
                    # pick some observables
                    for o, obs in enumerate(observable_list):
                        observable = extract_observable(case, obs)
                        observables[i,o] = observable
        # plot
        x_axis = get_x_axis(x_axis_list, x_axis_variable, case)
        inds = np.argsort(x_axis)
        x_axis = x_axis[inds]
        for o, obs in enumerate(observable_list):
            y_axis = observables[inds, o]
            if obs["axis"] == 1:
                ax1.plot(x_axis, y_axis, label = line + " " + obs["label"],
                    linestyle = obs["linestyle"], color = colors[j])
            if obs["axis"] == 2:
                ax2.plot(x_axis, y_axis, label = line + " " + obs["label"],
                    linestyle = obs["linestyle"], color = colors[j])
    ax1.legend()
    ax2.legend()
    plt.savefig(out_filename)        

def get_x_axis(x_axis_list, x_axis_variable, case):
    if x_axis_variable == "method":
        print("categorical data >:(")
    if x_axis_variable == "mesh":
        m = []
        for mesh_size in x_axis_list:
            if case["method"] == "meshed":
                match mesh_size:
                    case 2:
                        m.append(836)
                    case 3:
                        m.append(1246)
                    case 3:
                        m.append(1720)
            else:
                match mesh_size:
                    case 2:
                        m.append(1330)
                    case 3:
                        m.append(2110)
                    case 3:
                        m.append(2546)
        x_axis = np.array(m)
    if x_axis_variable == "Re":
        x_axis = np.array(x_axis_list)
    if x_axis_variable == "chi":
        x_axis = np.array(x_axis_list)
    if x_axis_variable == "radius":
        x_axis = np.array(x_axis_list)
    if x_axis_variable == "rmax":
        x_axis = np.array(x_axis_list)
    return x_axis

def extract_observable(case, observable):
    if observable["type"] == "force":
        result = extract_force_observable(case, observable["circle"], 
            observable["component"])
    return result

def extract_force_observable(case, circle, component):
    # this should probably be nans or something
    name = "force_" + circle
    data = case["observables"]
    result = 0.0
    for data_point in data:
        if data_point["name"] == name:
            if component in data_point["mean"]:
                result = data_point["mean"][component]
    return result
