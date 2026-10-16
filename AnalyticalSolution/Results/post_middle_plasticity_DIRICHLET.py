import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.colors as mcolors
import os

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",   # (Optional) you can force 'serif' font if you want classic LaTeX look
})

def adjust_color(rgb, amount=0.35, darker=True):
    """
    Schiarisce o scurisce un colore RGB di una percentuale `amount`.
    - darker=True  → avvicina il colore al nero   (più scuro)
    - darker=False → avvicina il colore al bianco (più chiaro)
    """
    c = np.array(mcolors.to_rgb(rgb))
    if darker:
        c = c * (1 - amount)          # verso il nero
    else:
        c = 1 - (1 - c) * (1 - amount)  # verso il bianco
    return tuple(np.clip(c, 0, 1))

def conv_plot(fig, ax, L2err, h, p, scheme, see_ref= False, label = "", L2_plot= False):
    """
    Disegna in scala log‑log l’andamento dell’errore L2 al variare di h
    e, opzionalmente, le rette di riferimento di ordine 1,2,…,n_degree.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Asse su cui tracciare il grafico.
    L2err : array‑like
        Valori dell’errore L2.
    h : array‑like
        Passi di mesh corrispondenti (stessa lunghezza di L2err).
    label : str, optional
        Etichetta per la curva dell’errore.
    n_degree : int, optional
        Numero di pendenze di riferimento da mostrare.
    see_ref : bool, optional
        Se True mostra le rette di riferimento.
    """
    # converto in array NumPy per sicurezza
    L2err = np.asarray(L2err, dtype=float)
    h      = np.asarray(h, dtype=float)

    # if (L2_plot):
    reference_p = 2
    # else:
    #     reference_p = p

    # rette di riferimento (ordine 1, 2, …)
    if see_ref:
        for k in range(reference_p, reference_p + 1):
            ref = (L2err[-1] * (h / h[-1])**k)*0.6         # stessa origine, pendenza k
            ax.plot(h, ref, "--", linewidth=linesize-0.2, color = 'black')

            tri_orig_linear = [h[-1], ref[-1]]
            draw_loglog_slope(fig, ax, tri_orig_linear, tri_width, k, inverted=False, polygon_kwargs=tri_kwargs, label_kwargs=tri_label_kwargs)

    # curva dell’errore numerico
    # ax.plot(h, L2err, color=color, linestyle='-', marker=marker, linewidth=linesize, markersize=markersize, alpha=0.8, label = label)
    st_p      = p_styles[p]
    st_scheme = scheme_styles[scheme]

    color = st_scheme["color"]
        
    ax.plot(h, L2err,
            color      = color,
            marker     = st_p["marker"],
            linestyle  = st_scheme["linestyle"],
            linewidth  = st_scheme["lw"],
            markersize = st_scheme["markersize"],
            alpha      = st_scheme["alpha"],
            label = label)


def draw_loglog_slope(fig, ax, origin, width_inches, slope, inverted=False, color=None, polygon_kwargs=None, label=True, labelcolor=None, label_kwargs=None, zorder=None):
    """
    This function draws slopes or "convergence triangles" into loglog plots.
    @param fig: The figure
    @param ax: The axes object to draw to
    @param origin: The 2D origin (usually lower-left corner) coordinate of the triangle
    @param width_inches: The width in inches of the triangle
    @param slope: The slope of the triangle, i.e. order of convergence
    @param inverted: Whether to mirror the triangle around the origin, i.e. whether
        it indicates the slope towards the lower left instead of upper right (defaults to false)
    @param color: The color of the of the triangle edges (defaults to default color)
    @param polygon_kwargs: Additional kwargs to the Polygon draw call that creates the slope
    @param label: Whether to enable labeling the slope (defaults to true)
    @param labelcolor: The color of the slope labels (defaults to the edge color)
    @param label_kwargs: Additional kwargs to the Annotation draw call that creates the labels
    @param zorder: The z-order value of the triangle and labels, defaults to a high value
    """

    if polygon_kwargs is None:
        polygon_kwargs = {}
    if label_kwargs is None:
        label_kwargs = {}

    if color is not None:
        polygon_kwargs["color"] = color
    if "linewidth" not in polygon_kwargs:
        polygon_kwargs["linewidth"] = 0.75 * matplotlib.rcParams["lines.linewidth"]
    if labelcolor is not None:
        label_kwargs["color"] = labelcolor
    if "color" not in label_kwargs:
        label_kwargs["color"] = polygon_kwargs["color"]
    if "fontsize" not in label_kwargs:
        label_kwargs["fontsize"] = 1.2 * matplotlib.rcParams["font.size"]

    if inverted:
        width_inches = -width_inches
    if zorder is None:
        zorder = 10

    # For more information on coordinate transformations in Matplotlib see
    # https://matplotlib.org/3.1.1/tutorials/advanced/transforms_tutorial.html

    # Convert the origin into figure coordinates in inches
    origin_disp = ax.transData.transform(origin)
    origin_dpi = fig.dpi_scale_trans.inverted().transform(origin_disp)

    # Obtain the bottom-right corner in data coordinates
    corner_dpi = origin_dpi + width_inches * np.array([1.0, 0.0])
    corner_disp = fig.dpi_scale_trans.transform(corner_dpi)
    corner = ax.transData.inverted().transform(corner_disp)

    (x1, y1) = (origin[0], origin[1])
    x2 = corner[0]

    # The width of the triangle in data coordinates
    width = x2 - x1
    # Compute offset of the slope
    log_offset = y1 / (x1 ** slope)

    y2 = log_offset * ((x1 + width) ** slope)
    height = y2 - y1

    # The vertices of the slope
    a = origin
    b = corner
    c = [x2, y2]

    # Draw the slope triangle
    X = np.array([a, b, c])
    triangle = plt.Polygon(X[:3,:], fill=False, zorder=zorder, **polygon_kwargs)
    ax.add_patch(triangle)

    # Convert vertices into display space
    a_disp = ax.transData.transform(a)
    b_disp = ax.transData.transform(b)
    c_disp = ax.transData.transform(c)

    # Figure out the center of the triangle sides in display space
    bottom_center_disp = a_disp + 0.5 * (b_disp - a_disp)
    bottom_center = ax.transData.inverted().transform(bottom_center_disp)

    right_center_disp = b_disp + 0.5 * (c_disp - b_disp)
    right_center = ax.transData.inverted().transform(right_center_disp)

    # Label alignment depending on inversion parameter
    va_xlabel = "top" if not inverted else "bottom"
    ha_ylabel = "left" if not inverted else "right"

    # Label offset depending on inversion parameter
    offset_xlabel = [0.0, -0.33 * label_kwargs["fontsize"]] if not inverted else [0.0, 0.33 * label_kwargs["fontsize"]]
    offset_ylabel = [0.33 * label_kwargs["fontsize"], 0.0] if not inverted else [-0.33 * label_kwargs["fontsize"], 0.0]

    # Draw the slope labels
    ax.annotate(r'$1$', bottom_center, xytext=offset_xlabel, textcoords='offset points', ha="center", va=va_xlabel, zorder=zorder, **label_kwargs)
    ax.annotate(fr'${slope}$', right_center, xytext=offset_ylabel, textcoords='offset points', ha=ha_ylabel, va="center", zorder=zorder, **label_kwargs)


plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
linesize = 1.2
markersize = 3
fontsize = 15
color1 = (0.2, 0.4, 1.0)
color2 = (0.8, 0.2, 0.2)
color4 = (0.2, 0.7, 0.2)  
color3 = (0.6, 0.2, 0.8) 
color5 = (1.0, 0.5, 0.0)   # Orange
color6 = (0.0, 0.6, 0.6)   # Teal (blue-green)

p_styles = {
    1: {"color": (0.20, 0.40, 1.00), "marker": "o"},   # blu,  cerchio
    2: {"color": (0.80, 0.20, 0.20), "marker": "s"},   # rosso, quadrato
    3: {"color": (0.20, 0.70, 0.20), "marker": "^"},   # verde, rombo
    4: {"color": (0.60, 0.20, 0.80), "marker": "^"},   # viola, triangolo ^
    5: {"color": (1.00, 0.50, 0.00), "marker": "v"}   # arancio, triangolo v
}

scheme_styles = {
    "gap_sbm":      {"linestyle": "-",  "alpha": 0.95, "lw": linesize, "markersize":markersize, "color": color1},
    "sbm": {"linestyle": "--", "alpha": 0.95, "lw": linesize, "markersize":markersize, "color": color1},
    "standard": {"linestyle": "--",  "alpha": 0.95, "lw": linesize, "markersize":markersize, "color": color2},
    "body-fitted": {"linestyle": "-", "alpha": 0.8, "lw": linesize*1.5, "markersize":markersize, "color": "#FFCC00"},
}


# # Draw convergence triangles
tri_width = 0.015
tri_kwargs = {
    "color" : 'k',
    "alpha" : 0.7
}
tri_label_kwargs = {
    "color" : 'k',
    "alpha" : 0.7
}


# Mixed derivatives in SBM operator
fig_L2, ax_L2 = plt.subplots(1, 1, figsize=(6, 7))  # Two subplots for norm L2 and H1
fig_H1, ax_H1 = plt.subplots(1, 1, figsize=(6, 7))

### p1
h_list = [0.8, 0.4, 0.2, 0.1, 0.05, 0.025]
L2_error = [0.2686022480387245, 0.0804386701568839, 0.023261680188645065, 0.004371007114721615, 0.0013611263125460414, 0.0003716415962316471]
H1_error = [0.3695021822911408, 0.13367250906512043, 0.053401959693192204, 0.02185659174313961, 0.010437788404208897, 0.005074628519952715]
conv_plot(fig_L2, ax_L2, L2_error, h_list, 1, scheme="standard", see_ref=True, L2_plot= True, label = r"body-fitted")
conv_plot(fig_H1, ax_H1, H1_error, h_list, 1, scheme="standard", see_ref=True, label = r"body-fitted")


h_list = [0.9400000000000001, 0.47000000000000003, 0.23500000000000001, 0.11750000000000001, 0.058750000000000004, 0.029375000000000002]
L2_error = [0.2183653850022006, 0.07134348698257013, 0.026371976308893083, 0.007022283983969592, 0.001952869732814172, 0.00046216225130706236]
H1_error = [0.3049457281567717, 0.14788587713978027, 0.06298262631059666, 0.027171274326510425, 0.01268036439879562, 0.006089084375163752]
conv_plot(fig_L2, ax_L2, L2_error, h_list, 1, scheme="sbm", see_ref=True, L2_plot= True, label = r"sbm")
conv_plot(fig_H1, ax_H1, H1_error, h_list, 1, scheme="sbm", see_ref=True, label = r"sbm")

# ### p2  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# body fitted
h_list = [0.8, 0.4, 0.2, 0.1, 0.05, 0.025]
L2_error = [0.015027702312164097, 0.0022270954124757446, 0.00038491290901642377, 0.0009935530858439766, 0.00017461312880171954, 4.613845931006117e-05]
H1_error = [0.029409461579165068, 0.007276167737756561, 0.0016184568873353897, 0.0020679809538110727, 0.00046706537653014296, 0.00016505280086486531]

conv_plot(fig_L2, ax_L2, L2_error, h_list, 2, scheme="standard", see_ref=True, L2_plot= True)
conv_plot(fig_H1, ax_H1, H1_error, h_list, 2, scheme="standard", see_ref=True)

#sbm

h_list = [0.9400000000000001, 0.47000000000000003, 0.23500000000000001, 0.11750000000000001, 0.058750000000000004, 0.029375000000000002]
L2_error = [0.024409722008257823, 0.006749922169950314, 0.0011638541884788779, 0.0006406705629395438, 9.249074537055336e-05, 6.516351392623212e-05]
H1_error = [0.03864385395183968, 0.009750294173159415, 0.002480014641727451, 0.0017326372002831206, 0.00043943313061730157, 0.00024115259809702754]
conv_plot(fig_L2, ax_L2, L2_error, h_list, 2, scheme="sbm", see_ref=False, L2_plot= True)
conv_plot(fig_H1, ax_H1, H1_error, h_list, 2, scheme="sbm", see_ref=False)


# # last
# tri_orig_linear = [h[-2], yDegreesLast[-1]]
# draw_loglog_slope(fig, axs, tri_orig_linear, tri_width, i, inverted=False, polygon_kwargs=tri_kwargs, label_kwargs=tri_label_kwargs)

ax_L2.set_xscale('log')
ax_L2.set_yscale('log')
ax_L2.grid(True, which='both', color='lightgrey', linewidth=0.8)

ax_L2.set_xlabel('$h$', fontsize=fontsize)
ax_L2.set_ylabel('$\\| {u}-{u_h}\\|_{L^2} / \\| {u}\\|_{L^2}$', fontsize=fontsize)

ax_H1.set_xscale('log')
ax_H1.set_yscale('log')
ax_H1.grid(True, which='both', color='lightgrey', linewidth=0.8)

ax_H1.set_xlabel('$h$', fontsize=fontsize)
ax_H1.set_ylabel('$\\| {u}-{u_h}\\|_{H^1} / \\| {u}\\|_{H^1}$', fontsize=fontsize)

from matplotlib.lines import Line2D

# --- handle per l’ordine p (come già facevi) ---
handles_p = [Line2D([], [], color=st["color"],
                    marker=st["marker"], linestyle='',
                    markersize=markersize*1.5,
                    label=f"$p={p}$")
             for p, st in p_styles.items()]

# --- handle per la retta di riferimento ---
handle_ref = Line2D([], [], color='black', linestyle='--',
                    linewidth=linesize-0.2, label='Reference')

# --- legenda unica ---
# all_handles = handles_p + handles_scheme + [handle_ref]
# all_handles = handles_p #+ [handle_ref]
# ax_L2.legend(handles=all_handles, title="", frameon=False,
#              loc="lower right", fontsize=fontsize)

ax_L2.set_title(r"SBM Vs GAP-SBM", fontsize=fontsize)
ax_L2.legend(loc='lower right', frameon=False, fontsize=fontsize*0.8)
ax_H1.set_title(r"SBM Vs GAP-SBM ($H^1$)", fontsize=fontsize)
ax_H1.legend(loc='lower right', frameon=False, fontsize=fontsize*0.8)

# Adjust layout and save plot with name matching this folder
fig_L2.tight_layout()
fig_H1.tight_layout()
# Derive PDF name from this script's folder (replace spaces with underscores)
_folder_name = os.path.basename(os.path.dirname(__file__))
_pdf_name_L2 = _folder_name.replace(' ', '_') + ".pdf"
_pdf_name_H1 = _folder_name.replace(' ', '_') + "_H1.pdf"
fig_L2.savefig(_pdf_name_L2, format="pdf", dpi=1000, bbox_inches="tight")
fig_H1.savefig(_pdf_name_H1, format="pdf", dpi=1000, bbox_inches="tight")

plt.show()
