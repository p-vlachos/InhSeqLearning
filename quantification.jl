# using PyPlot
# using ImageFiltering
# using LinearAlgebra
# using StatsBase
# using CurveFit


# Nmaxmembers = 300
# Ni2 = 500
# dt = .1
# Npop = size(popmembers)[1]
# Npop == 9 ? Nseq = 3 : Nseq = 5
# Nsteps = round(Int, T/dt)
# fr_threshold = 4   # Spiking threshold for assemblies (Hz; Chenkov:30 spikes/sec)
# fRate_window = 5
# time_limit = 200
# overlap_limit = 20
# Npop == 9 ? seq_length = 3 : seq_length = 4

function plot_eeWeights(popmembers, Npop)
    # Plot heatmap of the total average E-to-E assemblies synaptic coupling strength
    avg_weights = zeros(Npop, Npop)
    labels = String[]
    label_size = 13
    axis_size = 11
    for ipop = 1:Npop
        imembers = filter(i->(i>0), popmembers[ipop,:])
        append!(labels, [string(ipop)])
        for jpop = 1:Npop
            jmembers = filter(i->(i>0), popmembers[jpop, :])
            avg_weights[ipop, jpop] = sum(new_weights[imembers, jmembers]) / count(i->(i>0), new_weights[imembers, jmembers])
        end
    end
    top_cmap = matplotlib.cm.get_cmap("Blues", 128)
    bottom_cmap = matplotlib.cm.get_cmap("gray", 128)
    newcolors = matplotlib.colors.ListedColormap([bottom_cmap(1:128); top_cmap(1:128)], name="BlGy")
    matshow(avg_weights, origin="lower", cmap=newcolors)
    xlabel("Postsynaptic E-assembly", fontsize=label_size)
    ylabel("Presynaptic E-assembly", fontsize=label_size)
    cb = colorbar(shrink=0.8, pad=.1)
    cb_ticks = matplotlib.ticker.MaxNLocator(nbins=5); cb.locator = cb_ticks; cb.update_ticks()
    cb.set_label(label="Coupling Strength (pF)", size=label_size)
    cb.ax.tick_params(labelsize=axis_size)
    ax = gca()
    major_ticks = matplotlib.ticker.FixedLocator([1.5, 5.5, 9.5, 13.5, 17.5])
    minor_ticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_ticks); ax.yaxis.set_major_locator(major_ticks)
    ax.xaxis.set_minor_locator(minor_ticks); ax.yaxis.set_minor_locator(minor_ticks)
    ax.set_xticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray", fontfamily="monospace")
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.set_yticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray", fontfamily="monospace")
    ax.set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=true, labelleft=false, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    ax.axhline([3.5], 0, .4, linestyle="-", linewidth=.75, color="tab:red"); ax.axhline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75, color="tab:red");
    ax.axhline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75, color="tab:red"); ax.axhline([15.5], 0.6, 1, linestyle="-", linewidth=.75, color="tab:red");
    ax.axhline([3.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axhline([7.5], linestyle=":", linewidth=.75, color="tab:red");
    ax.axhline([11.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axhline([15.5], linestyle=":", linewidth=.75, color="tab:red");
    ax.axvline([3.5], 0, .4, linestyle="-", linewidth=.75, color="tab:red"); ax.axvline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75, color="tab:red");
    ax.axvline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75, color="tab:red"); ax.axvline([15.5], 0.6, 1, linestyle="-", linewidth=.75, color="tab:red");
    ax.axvline([3.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axvline([7.5], linestyle=":", linewidth=.75, color="tab:red");
    ax.axvline([11.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axvline([15.5], linestyle=":", linewidth=.75, color="tab:red");
    text(4, 23, "Sequence Label", fontsize=label_size)
    tight_layout()
    savefig("avgE-to-E_strength.png", bbox_inches="tight", dpi=150)
    PyPlot.clf()
end

function plot_eeWeights(new_weights, popmembers, Npop, path)
    # Plot heatmap of the total average E-to-E assemblies synaptic coupling strength
    avg_weights = zeros(Npop, Npop)
    labels = String[]
    label_size = 13
    axis_size = 11
    for ipop = 1:Npop
        imembers = filter(i->(i>0), popmembers[ipop,:])
        append!(labels, [string(ipop)])
        for jpop = 1:Npop
            jmembers = filter(i->(i>0), popmembers[jpop, :])
            avg_weights[ipop, jpop] = sum(new_weights[imembers, jmembers]) / count(i->(i>0), new_weights[imembers, jmembers])
        end
    end
    top_cmap = matplotlib.cm.get_cmap("Blues", 128)
    bottom_cmap = matplotlib.cm.get_cmap("gray", 128)
    newcolors = matplotlib.colors.ListedColormap([bottom_cmap(1:128); top_cmap(1:128)], name="BlGy")
    matshow(avg_weights, origin="lower", cmap=newcolors)
    xlabel("Postsynaptic E-assembly", fontsize=label_size)
    ylabel("Presynaptic E-assembly", fontsize=label_size)
    cb = colorbar(shrink=0.8, pad=.1)
    cb_ticks = matplotlib.ticker.MaxNLocator(nbins=5); cb.locator = cb_ticks; cb.update_ticks()
    cb.set_label(label="Coupling Strength (pF)", size=label_size)
    cb.ax.tick_params(labelsize=axis_size)
    ax = gca()
    major_ticks = matplotlib.ticker.FixedLocator([1.5, 5.5, 9.5, 13.5, 17.5])
    minor_ticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_ticks); ax.yaxis.set_major_locator(major_ticks)
    ax.xaxis.set_minor_locator(minor_ticks); ax.yaxis.set_minor_locator(minor_ticks)
    ax.set_xticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray", fontfamily="monospace")
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.set_yticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray", fontfamily="monospace")
    ax.set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=true, labelleft=false, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    ax.axhline([3.5], 0, .4, linestyle="-", linewidth=.75, color="tab:red"); ax.axhline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75, color="tab:red");
    ax.axhline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75, color="tab:red"); ax.axhline([15.5], 0.6, 1, linestyle="-", linewidth=.75, color="tab:red");
    ax.axhline([3.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axhline([7.5], linestyle=":", linewidth=.75, color="tab:red");
    ax.axhline([11.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axhline([15.5], linestyle=":", linewidth=.75, color="tab:red");
    ax.axvline([3.5], 0, .4, linestyle="-", linewidth=.75, color="tab:red"); ax.axvline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75, color="tab:red");
    ax.axvline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75, color="tab:red"); ax.axvline([15.5], 0.6, 1, linestyle="-", linewidth=.75, color="tab:red");
    ax.axvline([3.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axvline([7.5], linestyle=":", linewidth=.75, color="tab:red");
    ax.axvline([11.5], linestyle=":", linewidth=.75, color="tab:red"); ax.axvline([15.5], linestyle=":", linewidth=.75, color="tab:red");
    text(4, 23, "Sequence Label", fontsize=label_size)
    tight_layout()
    savefig(string(path, "avgE-to-E_strength.png"), bbox_inches="tight", dpi=150)
    PyPlot.clf()
end

function plot_eiWeights(popmembers, ipopmembers, Npop)
    # Plot heatmap of the total average I₂-to-E assemblies synaptic coupling strength
    avg_weights = zeros(Npop, Npop)
    label_size = 13
    axis_size = 11
    for iipop = 1:Npop
        imembers = ipopmembers[iipop,:]
        for ipop = 1:Npop
            members = filter(i->(i>0), popmembers[ipop, :])
            avg_weights[iipop, ipop] = sum(new_weights[imembers, members]) / count(i->(i>0), new_weights[imembers, members])
        end
    end
    matshow(avg_weights, origin="lower", cmap="RdGy_r")
    xlabel("Postsynaptic E-assembly", fontsize=label_size)
    ylabel("Presynaptic I₂-assembly", fontsize=label_size)
    cb = colorbar(shrink=0.8, pad=.1)
    cb_ticks = matplotlib.ticker.MaxNLocator(nbins=5); cb.locator = cb_ticks; cb.update_ticks()
    cb.set_label(label="Coupling Strength (pF)", size=label_size)
    cb.ax.tick_params(labelsize=axis_size)
    ax = gca()
    major_ticks = matplotlib.ticker.FixedLocator([1.5, 5.5, 9.5, 13.5, 17.5])
    minor_ticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_ticks); ax.yaxis.set_major_locator(major_ticks)
    ax.xaxis.set_minor_locator(minor_ticks); ax.yaxis.set_minor_locator(minor_ticks)
    ax.set_xticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.set_yticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=true, labelleft=false, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    ax.axhline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axhline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axhline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axhline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axhline([3.5], linestyle=":", linewidth=.75); ax.axhline([7.5], linestyle=":", linewidth=.75);
    ax.axhline([11.5], linestyle=":", linewidth=.75); ax.axhline([15.5], linestyle=":", linewidth=.75);
    ax.axvline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axvline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axvline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axvline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axvline([3.5], linestyle=":", linewidth=.75); ax.axvline([7.5], linestyle=":", linewidth=.75);
    ax.axvline([11.5], linestyle=":", linewidth=.75); ax.axvline([15.5], linestyle=":", linewidth=.75);
    text(4, 23, "Sequence Label", fontsize=label_size)
    tight_layout()
    savefig("avgI2-to-E_strength.png", bbox_inches="tight", dpi=150)
    PyPlot.clf()
end

function plot_eiWeights(new_weights, popmembers, ipopmembers, Npop, path)
    # Plot heatmap of the total average I₂-to-E assemblies synaptic coupling strength
    avg_weights = zeros(Npop, Npop)
    label_size = 13
    axis_size = 11
    for iipop = 1:Npop
        imembers = ipopmembers[iipop,:]
        for ipop = 1:Npop
            members = filter(i->(i>0), popmembers[ipop, :])
            avg_weights[iipop, ipop] = sum(new_weights[imembers, members]) / count(i->(i>0), new_weights[imembers, members])
        end
    end
    matshow(avg_weights, origin="lower", cmap="RdGy_r")
    xlabel("Postsynaptic E-assembly", fontsize=label_size)
    ylabel("Presynaptic I₂-assembly", fontsize=label_size)
    cb = colorbar(shrink=0.8, pad=.1)
    cb_ticks = matplotlib.ticker.MaxNLocator(nbins=5); cb.locator = cb_ticks; cb.update_ticks()
    cb.set_label(label="Coupling Strength (pF)", size=label_size)
    cb.ax.tick_params(labelsize=axis_size)
    ax = gca()
    major_ticks = matplotlib.ticker.FixedLocator([1.5, 5.5, 9.5, 13.5, 17.5])
    minor_ticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_ticks); ax.yaxis.set_major_locator(major_ticks)
    ax.xaxis.set_minor_locator(minor_ticks); ax.yaxis.set_minor_locator(minor_ticks)
    ax.set_xticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.set_yticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=true, labelleft=false, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    ax.axhline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axhline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axhline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axhline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axhline([3.5], linestyle=":", linewidth=.75); ax.axhline([7.5], linestyle=":", linewidth=.75);
    ax.axhline([11.5], linestyle=":", linewidth=.75); ax.axhline([15.5], linestyle=":", linewidth=.75);
    ax.axvline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axvline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axvline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axvline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axvline([3.5], linestyle=":", linewidth=.75); ax.axvline([7.5], linestyle=":", linewidth=.75);
    ax.axvline([11.5], linestyle=":", linewidth=.75); ax.axvline([15.5], linestyle=":", linewidth=.75);
    text(4, 23, "Sequence Label", fontsize=label_size)
    tight_layout()
    savefig(string(path, "avgI2-to-E_strength.png"), bbox_inches="tight", dpi=150)
    PyPlot.clf()
end

function plot_eiWeights_zoom(popmembers, ipopmembers, Npop; seq_num=1, seq_len=4)
    # Plot heatmap of the total average I₂-to-E assemblies synaptic coupling strength
    # (zoomed in to the pre-synaptic assembly 'seq_num')
    endInd = seq_num*seq_len
    startInd = endInd - (seq_len-1)
    avg_weights = zeros(seq_len, Npop)  # Contains the average strength from I-to-E
    label_size = 16
    axis_size = 14
    line_width = 3.
    for iipop = startInd:endInd
        (mod(iipop, seq_len) == 0) ? index = seq_len : index = mod(iipop, seq_len)
        imembers = ipopmembers[iipop,:]
        for ipop = 1:Npop
            members = filter(i->(i>0), popmembers[ipop, :])
            avg_weights[index, ipop] = sum(new_weights[imembers, members]) / count(i->(i>0), new_weights[imembers, members])
        end
    end
    matshow(avg_weights, origin="lower", cmap="RdGy_r")
    xlabel("Postsynaptic E-assembly", fontsize=label_size)
    ylabel("Presynaptic\nI₂-assembly", fontsize=label_size)
    cb = colorbar(shrink=0.8, pad=.01)
    cb_ticks = matplotlib.ticker.MaxNLocator(nbins=5); cb.locator = cb_ticks; cb.update_ticks()
    cb.set_label(label="Coupling\nStrength (pF)", size=label_size)
    cb.ax.tick_params(labelsize=axis_size)
    ax = gca()
    major_xticks = matplotlib.ticker.FixedLocator([1.5, 5.5, 9.5, 13.5, 17.5])
    minor_xticks = matplotlib.ticker.MultipleLocator(1)
    major_yticks = matplotlib.ticker.FixedLocator([1.5])
    minor_yticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_xticks); ax.yaxis.set_major_locator(major_yticks)
    ax.xaxis.set_minor_locator(minor_xticks); ax.yaxis.set_minor_locator(minor_yticks)
    ax.set_xticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=label_size, color="tab:gray")
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.set_yticklabels(["", startInd, "", ""], minor=true, fontsize=axis_size, color="tab:gray")
    if seq_num == 1
        ax.set_yticklabels(["A"], minor=false, fontsize=label_size, fontfamily="monospace")
        ax.axvline([3.5], linestyle="-", linewidth=line_width);
    elseif seq_num == 2
        ax.set_yticklabels(["B"], minor=false, fontsize=label_size, fontfamily="monospace")
        ax.axvline([3.5], linestyle="-", linewidth=line_width);ax.axvline([7.5], linestyle="-", linewidth=line_width)
    elseif seq_num == 3
        ax.set_yticklabels(["C"], minor=false, fontsize=label_size, fontfamily="monospace")
        ax.axvline([7.5], linestyle="-", linewidth=line_width);ax.axvline([11.5], linestyle="-", linewidth=line_width)
    elseif seq_num == 4
        ax.set_yticklabels(["D"], minor=false, fontsize=label_size, fontfamily="monospace")
        ax.axvline([11.5], linestyle="-", linewidth=line_width);ax.axvline([15.5], linestyle="-", linewidth=line_width);
    elseif seq_num == 5
        ax.set_yticklabels(["E"], minor=false, fontsize=label_size, fontfamily="monospace")
        ax.axvline([15.5], linestyle="-", linewidth=line_width);
    end
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    ax.axvline([3.5], linestyle=":", linewidth=line_width); ax.axvline([7.5], linestyle=":", linewidth=line_width);
    ax.axvline([11.5], linestyle=":", linewidth=line_width); ax.axvline([15.5], linestyle=":", linewidth=line_width);
    text(7.9, 4.5, "Sequence Label", fontsize=label_size)
    tight_layout()
    savefig("avgI2-to-E_strength_zoom.png", dpi=150, bbox_inches="tight")
    PyPlot.clf()
end

function plot_eiWeights_zoom2(popmembers, ipopmembers, Npop; seq_num=1, seq_len=4)
    # Plot heatmap of the total average I₂-to-E assemblies synaptic coupling strength
    # (zoomed-in to E/I-assembly 'seq_num')
    endInd = seq_num*seq_len
    startInd = endInd - (seq_len-1)
    avg_weights = zeros(seq_len, seq_len)
    for iipop = startInd:endInd
        (mod(iipop, seq_len) == 0) ? index = seq_len : index = mod(iipop, seq_len)
        imembers = ipopmembers[iipop,:]
        for ipop = startInd:endInd
            (mod(ipop, seq_len) == 0) ? indexi = seq_len : indexi = mod(ipop, seq_len)
            members = filter(i->(i>0), popmembers[ipop, :])
            avg_weights[index, indexi] = sum(new_weights[imembers, members]) / count(i->(i>0), new_weights[imembers, members])
        end
    end
    matshow(avg_weights, origin="lower", cmap="RdGy_r")
    xlabel("Postsynaptic E-assembly")
    ylabel("Presynaptic I₂-assembly")
    # colorbar(shrink=0.8, label="Coupling Strength (pF)")
    ax = gca()
    major_ticks = matplotlib.ticker.FixedLocator([1.5])
    minor_ticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_ticks); ax.yaxis.set_major_locator(major_ticks)
    ax.xaxis.set_minor_locator(minor_ticks); ax.yaxis.set_minor_locator(minor_ticks)
    ax.set_xticklabels(["", startInd, "", ""], minor=true, fontsize="x-small", color="tab:gray")
    ax.set_yticklabels(["", startInd, "", ""], minor=true, fontsize="x-small", color="tab:gray")
    if seq_num == 1
        ax.set_yticklabels(["A"], minor=false, fontsize="large", fontfamily="monospace")
        ax.set_xticklabels(["A"], minor=false, fontsize="large", fontfamily="monospace")
    elseif seq_num == 2
        ax.set_yticklabels(["B"], minor=false, fontsize="large", fontfamily="monospace")
        ax.set_xticklabels(["B"], minor=false, fontsize="large", fontfamily="monospace")
    elseif seq_num == 3
        ax.set_yticklabels(["C"], minor=false, fontsize="large", fontfamily="monospace")
        ax.set_xticklabels(["C"], minor=false, fontsize="large", fontfamily="monospace")
    elseif seq_num == 4
        ax.set_yticklabels(["D"], minor=false, fontsize="large", fontfamily="monospace")
        ax.set_xticklabels(["D"], minor=false, fontsize="large", fontfamily="monospace")
    elseif seq_num == 5
        ax.set_yticklabels(["E"], minor=false, fontsize="large", fontfamily="monospace")
        ax.set_xticklabels(["E"], minor=false, fontsize="large", fontfamily="monospace")
    end
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    for (ind, val) in pairs(avg_weights)
        ax.text(ind[2]-1, ind[1]-1, round(val, digits=2), ha="center", va="center", fontsize="small", weight="medium", color="black")
        ax.add_patch(matplotlib.patches.Rectangle((ind[1]-1-0.25, ind[2]-1-0.06), 0.49, 0.14, color="white", fill="true", alpha=0.75))
    end
    tight_layout()
    savefig("avgI2-to-E_strength_zoom2.png")
    PyPlot.clf()
end

function plot_ieWeights(popmembers, ipopmembers, Npop)
    # Plot heatmap of weights from 2nd i-population assemblies to excitatory ones
    avg_weights = zeros(Npop, Npop)  # Contains the average strength from E-to-I2
    label_size = 13
    axis_size = 11
    for ipop = 1:Npop
        members = filter(i->(i>0), popmembers[ipop,:])
        for iipop = 1:Npop
            imembers = ipopmembers[iipop, :]
            avg_weights[ipop, iipop] = sum(new_weights[members, imembers]) / count(i->(i>0), new_weights[members, imembers])
        end
    end
    matshow(avg_weights, origin="lower", cmap="RdGy_r")
    xlabel("Postsynaptic I₂-assembly", size=label_size)
    ylabel("Presynaptic E-assembly", size=label_size)
    cb = colorbar(shrink=0.8, pad=.1)
    cb_ticks = matplotlib.ticker.MaxNLocator(nbins=5); cb.locator = cb_ticks; cb.update_ticks()
    cb.set_label(label="Coupling Strength (pF)", size=label_size)
    cb.ax.tick_params(labelsize=axis_size)
    ax = gca()
    major_ticks = matplotlib.ticker.FixedLocator([1.5, 5.5, 9.5, 13.5, 17.5])
    minor_ticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_ticks); ax.yaxis.set_major_locator(major_ticks)
    ax.xaxis.set_minor_locator(minor_ticks); ax.yaxis.set_minor_locator(minor_ticks)
    ax.set_xticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.set_yticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=true, labelleft=false, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    ax.axhline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axhline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axhline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axhline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axhline([3.5], linestyle=":", linewidth=.75); ax.axhline([7.5], linestyle=":", linewidth=.75);
    ax.axhline([11.5], linestyle=":", linewidth=.75); ax.axhline([15.5], linestyle=":", linewidth=.75);
    ax.axvline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axvline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axvline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axvline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axvline([3.5], linestyle=":", linewidth=.75); ax.axvline([7.5], linestyle=":", linewidth=.75);
    ax.axvline([11.5], linestyle=":", linewidth=.75); ax.axvline([15.5], linestyle=":", linewidth=.75);
    text(4, 23, "Sequence Label", fontsize=label_size)
    savefig("avgE-to-I2_strength.png", dpi=150, bbox_inches="tight")
    PyPlot.clf()
end

function plot_ieWeights(new_weights, popmembers, ipopmembers, Npop, path)
    # Plot heatmap of weights from 2nd i-population assemblies to excitatory ones
    avg_weights = zeros(Npop, Npop)  # Contains the average strength from E-to-I2
    label_size = 13
    axis_size = 11
    for ipop = 1:Npop
        members = filter(i->(i>0), popmembers[ipop,:])
        for iipop = 1:Npop
            imembers = ipopmembers[iipop, :]
            avg_weights[ipop, iipop] = sum(new_weights[members, imembers]) / count(i->(i>0), new_weights[members, imembers])
        end
    end
    matshow(avg_weights, origin="lower", cmap="RdGy_r")
    xlabel("Postsynaptic I₂-assembly", size=label_size)
    ylabel("Presynaptic E-assembly", size=label_size)
    cb = colorbar(shrink=0.8, pad=.1)
    cb_ticks = matplotlib.ticker.MaxNLocator(nbins=5); cb.locator = cb_ticks; cb.update_ticks()
    cb.set_label(label="Coupling Strength (pF)", size=label_size)
    cb.ax.tick_params(labelsize=axis_size)
    ax = gca()
    major_ticks = matplotlib.ticker.FixedLocator([1.5, 5.5, 9.5, 13.5, 17.5])
    minor_ticks = matplotlib.ticker.MultipleLocator(1)
    ax.xaxis.set_major_locator(major_ticks); ax.yaxis.set_major_locator(major_ticks)
    ax.xaxis.set_minor_locator(minor_ticks); ax.yaxis.set_minor_locator(minor_ticks)
    ax.set_xticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.set_yticklabels(["", 1, "", "", "", 5, "", "", "", 9, "", "", "", 13, "", "", "", 17], minor=true, fontsize=axis_size, color="tab:gray")
    ax.set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axis_size, fontfamily="monospace")
    ax.xaxis.set_tick_params(which="major", labeltop=true, labelbottom=false, bottom=false, top=false)
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, bottom=true, top=false)
    ax.yaxis.set_tick_params(which="major", labelright=true, labelleft=false, right=false, left=false)
    ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=true, right=false, left=true)
    ax.axhline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axhline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axhline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axhline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axhline([3.5], linestyle=":", linewidth=.75); ax.axhline([7.5], linestyle=":", linewidth=.75);
    ax.axhline([11.5], linestyle=":", linewidth=.75); ax.axhline([15.5], linestyle=":", linewidth=.75);
    ax.axvline([3.5], 0, .4, linestyle="-", linewidth=.75); ax.axvline([7.5], 0.2, 0.6, linestyle="-", linewidth=.75);
    ax.axvline([11.5], 0.4, 0.8, linestyle="-", linewidth=.75); ax.axvline([15.5], 0.6, 1, linestyle="-", linewidth=.75);
    ax.axvline([3.5], linestyle=":", linewidth=.75); ax.axvline([7.5], linestyle=":", linewidth=.75);
    ax.axvline([11.5], linestyle=":", linewidth=.75); ax.axvline([15.5], linestyle=":", linewidth=.75);
    text(4, 23, "Sequence Label", fontsize=label_size)
    savefig(string(path, "avgE-to-I2_strength.png"), dpi=150, bbox_inches="tight")
    PyPlot.clf()
end

function barplot_eiRecWeights(popmembers, ipopmembers, Npop)
    avg_weights = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers[pp, :])
        imembers = ipopmembers[pp, :]
        avg_weights[pp] = sum(new_weights[imembers, emembers]) / count(i->i>0,new_weights[imembers, emembers])
    end
    x = 1:Npop
    figure(figsize=(8, 4))
    xlim(0, Npop+1)
    bar(x, avg_weights, color="tab:red", width=0.5)
    ax = gca()
    ax.hlines(mean(avg_weights), 0, (Npop+1), label="average", ls="--", color="tab:gray")
    ax.vlines([4.5, 8.5, 12.5, 16.5], 0, 100, label="average", ls=":", color="tab:blue")
    xlabel("E\\I₂ - Assembly Index")
    ylabel("Coupling Strength (pF)")
    xticks = matplotlib.ticker.FixedLocator(1:seq_length:(Npop-seq_length+1))
    ax.xaxis.set_minor_locator(xticks)
    ax.set_xticklabels(["1", "5", "9", "13", "17"], minor=true, fontsize="small", c="tab:gray")
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, top=false, bottom=true)
    xticks = matplotlib.ticker.FixedLocator((seq_length/2)+0.5:seq_length:Npop)
    ax.xaxis.set_major_locator(xticks)
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false)
    savefig("avgI2-E_strength.png", bbox_inches="tight", dpi=150)
    PyPlot.clf()
end

function barplot_ieRecWeights(popmembers, ipopmembers, Npop)
    avg_weights = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers[pp, :])
        imembers = ipopmembers[pp, :]
        avg_weights[pp] = sum(new_weights[emembers, imembers]) / count(i->i>0,new_weights[emembers, imembers])
    end
    x = 1:Npop
    figure(figsize=(8, 4))
    xlim(0, Npop+1)
    bar(x, avg_weights, color="tab:blue", width=0.5)
    ax = gca()
    ax.hlines(mean(avg_weights), 0, (Npop+1), label="average", ls="--", color="tab:gray")
    ax.vlines([4.5, 8.5, 12.5, 16.5], 0, 1.5, label="average", ls=":", color="tab:red")
    xlabel("E\\I₂ - Assembly Index")
    ylabel("Coupling Strength (pF)")
    xticks = matplotlib.ticker.FixedLocator(1:seq_length:(Npop-seq_length+1))
    ax.xaxis.set_minor_locator(xticks)
    ax.set_xticklabels(["1", "5", "9", "13", "17"], minor=true, fontsize="small", c="tab:gray")
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, top=false, bottom=true)
    xticks = matplotlib.ticker.FixedLocator((seq_length/2)+0.5:seq_length:Npop)
    ax.xaxis.set_major_locator(xticks)
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false)
    savefig("avgE-I2_strength.png")
    PyPlot.clf()
end

function compare_eiRecWeights(popmembers, ipopmembers, Npop, popmembers2, ipopmembers2)
    # THIS IS FOR E-to-I2 weights
    avg_weights = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers[pp, :])
        imembers = ipopmembers[pp, :]
        avg_weights[pp] = sum(new_weights[emembers, imembers]) / count(i->i>0,new_weights[emembers, imembers])
    end
    avg_weights2 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers2[pp, :])
        imembers = ipopmembers2[pp, :]
        avg_weights2[pp] = sum(w8si[emembers, imembers]) / count(i->i>0,w8si[emembers, imembers])
    end
    x = 1:Npop
    figure(figsize=(8, 4))
    xlim(0, Npop+1)
    ylim(0, 4)
    bar(x, avg_weights2, color="tab:cyan", width=0.5)
    bar((x.+0.2), avg_weights, color="tab:blue", width=0.5)
    ax = gca()
    ax.legend(["without I₁→E plasticity", "with I₁→E plasticity"], fontsize="medium", loc="upper right", ncol=2)
    ax.hlines(mean(avg_weights), 0, (Npop+1), label="average", ls="--", color="tab:blue")
    ax.hlines(mean(avg_weights2), 0, (Npop+1), label="average", ls="--", color="tab:cyan")
    ax.vlines([4.6, 8.6, 12.6, 16.6], 0, 1.2, label="average", ls=":", color="tab:gray")
    xlabel("E\\I₂ - Assembly Index")
    ylabel("Average E→I₂ Coupling Strength (pF)")
    xticks = matplotlib.ticker.FixedLocator(1:seq_length:(Npop-seq_length+1))
    ax.xaxis.set_minor_locator(xticks)
    ax.set_xticklabels(["1", "5", "9", "13", "17"], minor=true, fontsize="small", c="tab:gray")
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, top=false, bottom=true)
    xticks = matplotlib.ticker.FixedLocator((seq_length/2)+0.5:seq_length:Npop)
    ax.xaxis.set_major_locator(xticks)
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false)
    savefig("avgE-I2_strength.png")
    PyPlot.clf()
end

function compare_ieRecWeights(popmembers, ipopmembers, Npop, popmembers2, ipopmembers2)
    # THIS IS FOR I2-to-E weights
    avg_weights = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers[pp, :])
        imembers = ipopmembers[pp, :]
        avg_weights[pp] = sum(new_weights[imembers, emembers]) / count(i->i>0,new_weights[imembers, emembers])
    end
    avg_weights2 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers2[pp, :])
        imembers = ipopmembers2[pp, :]
        avg_weights2[pp] = sum(w8si[imembers, emembers]) / count(i->i>0,w8si[imembers, emembers])
    end
    x = 1:Npop
    figure(figsize=(8, 4))
    xlim(0, Npop+1)
    ylim(0, 360)
    bar(x, avg_weights2, color="orangered", width=0.5)
    bar((x.+0.2), avg_weights, color="brown", width=0.5)
    ax = gca()
    # ax.legend(["without I₁→E plasticity", "with I₁→E plasticity"], fontsize="medium", loc="upper right", ncol=2)
    ax.legend(["without I₁", "with I₁"], fontsize="large", loc="upper right", ncol=2)
    ax.hlines(mean(avg_weights), 0, (Npop+1), label="average", ls="--", color="brown")
    ax.hlines(mean(avg_weights2), 0, (Npop+1), label="average", ls="--", color="orangered")
    ax.vlines([4.6, 8.6, 12.6, 16.6], 0, 120, label="average", ls=":", color="tab:gray")
    xlabel("E\\I₂ - Assembly Index", fontsize="x-large")
    ylabel("Average I₂→E Coupling Strength (pF)", fontsize="x-large")
    xticks = matplotlib.ticker.FixedLocator(1:seq_length:(Npop-seq_length+1))
    ax.xaxis.set_minor_locator(xticks)
    ax.set_xticklabels(["1", "5", "9", "13", "17"], minor=true, fontsize="medium", c="tab:gray")
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, top=false, bottom=true)
    xticks = matplotlib.ticker.FixedLocator((seq_length/2)+0.5:seq_length:Npop)
    ax.xaxis.set_major_locator(xticks)
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false)
    yticks(fontsize="large")
    savefig("avgI2-E_strength.png", bbox_inches="tight")
    PyPlot.clf()
end

function compare_eRecWeights(popmembers, popmembers2, Npop)
    avg_weights = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers[pp, :])
        avg_weights[pp] = sum(new_weights[emembers, emembers]) / count(i->i>0,new_weights[emembers, emembers])
    end
    avg_weights2 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers2[pp, :])
        avg_weights2[pp] = sum(w8si[emembers, emembers]) / count(i->i>0,w8si[emembers, emembers])
    end
    x = 1:Npop
    figure(figsize=(8, 4))
    xlim(0, Npop+1)
    ylim(0, 17.4)
    bar(x, avg_weights, color="tab:blue", width=0.5)
    bar((x.+0.2), avg_weights2, color="tab:cyan", width=0.5)
    ax = gca()
    # ax.legend(["with I₁→E plasticity", "without I₁→E plasticity"], fontsize="medium", loc="upper right", ncol=2)
    ax.legend(["with I₁", "without I₁"], fontsize="large", loc="upper right", ncol=2)
    ax.hlines(mean(avg_weights), 0, (Npop+1), label="average", ls="--", color="tab:blue")
    ax.hlines(mean(avg_weights2), 0, (Npop+1), label="average", ls="--", color="tab:cyan")
    ax.vlines([4.6, 8.6, 12.6, 16.6], 0, 5, label="average", ls=":", color="tab:gray")
    xlabel("E - Assembly Index", fontsize="x-large")
    ylabel("Average E Coupling Strength (pF)", fontsize="x-large")
    xticks = matplotlib.ticker.FixedLocator(1:seq_length:(Npop-seq_length+1))
    ax.xaxis.set_minor_locator(xticks)
    ax.set_xticklabels(["1", "5", "9", "13", "17"], minor=true, fontsize="medium", c="tab:gray")
    ax.xaxis.set_tick_params(which="minor", labeltop=false, labelbottom=true, top=false, bottom=true)
    xticks = matplotlib.ticker.FixedLocator((seq_length/2)+0.5:seq_length:Npop)
    ax.xaxis.set_major_locator(xticks)
    ax.set_xticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false)
    yticks(fontsize="large")
    savefig("avgE-E_strength.png", bbox_inches="tight")
    PyPlot.clf()
end

function getAssemblyfRates(Npop, T, eAssemSpikes, iAssemSpikes)
    # Returns the total firing rate for each E and I assembly through time (ms)
    Nsteps = round(Int, T/dt)
    eAssemfRates = zeros(Npop, T)
    iAssemfRates = zeros(Npop, T)
    ii = 1
    for tt = 1:Nsteps
        (tt > 10) && (ii = Int(trunc(tt*dt)))
        for ipop = 1:Npop
            eAssemfRates[ipop, ii] += eAssemSpikes[ipop, tt]
            iAssemfRates[ipop, ii] += iAssemSpikes[ipop, tt]
        end
    end
    return eAssemfRates, iAssemfRates
end

function calc_overlap(popmembers, Ne, Npop; Nmaxmembers=300)
    # Calculate percentage of neurons that bellong to more than one E-assembly
    stats = zeros(Int, Ne)
    for i = 1:Npop
        for j = 1:Nmaxmembers
            ind = popmembers[i, j]
            (ind == 0 || ind == -1) && continue
            stats[ind] += 1
        end
    end
    return round((count(i->(i>1), stats) / Ne) * 100, digits=2)
end

function calc_assOverlap(popmembers, Ne, Npop)
    stats = zeros(Npop, Ne)
    for pp = 1:Npop
        members = popmembers[pp, :][popmembers[pp,:] .> 0]
        stats[pp, members] .= 1
    end
    overlap = zeros(Npop, Npop)
    for pp1 = 1:Npop
        members = popmembers[pp1, :][popmembers[pp1,:] .> 0]
        for pp2 = 1:Npop
            members2 = popmembers[pp2, :][popmembers[pp2,:] .> 0]
            overlap[pp1, pp2] = round(sum(stats[pp1, :] .* stats[pp2, :]) ./ (length(members)) .* 100, digits=2)
        end
    end
    overlap[overlap .== 100] .= 0.
    return overlap
end

function plot_AssemfRate(Npop, T, eAssemfRates, iAssemfRates; fRate_window=2, seq_len=4)
    # ___________ Plot Assembly firing rates ___________
    # NOTE: GAUSSIAN FILTER
    # gaus_kernel = ImageFiltering.Kernel.gaussian((fRate_window,))
    # x= 1:T
    # # Print excitatory assembly firing rates through time
    # plot_data = zeros(Npop, T)
    # for pp = 1:Npop
    #     plot_data[pp,:] = imfilter(eAssemfRates[pp,:], gaus_kernel)
    # end

    # NOTE: ALPHA FUNCTION
    x= 1000:6000
    plot_data = convolveSpikes(times, interval=x)

    labels = String[]
    xlim(0, T)
    ylim(0, maximum(plot_data))
    for pp = 1:Npop
        if mod(pp, seq_len) == 1 && pp > 1
            legend(labels)
            xlabel("Simulation time (ms)")
            ylabel("Assembly firing rate (Hz)")
            savefig(string("E-AssemfRate_", pp-seq_len, "-", pp-1, ".png"))
            PyPlot.clf()
            xlim(0, T)
            ylim(0, maximum(plot_data))
            labels = String[]
        end
        plot(x, plot_data[pp,:])
        append!(labels, [string("E-", pp)])
    end
    legend(labels)
    xlabel("Simulation time (ms)")
    ylabel("E-Assembly firing rate (Hz)")
    savefig(string("E-AssemfRate_", Npop-seq_len, "-", Npop, ".png"))
    PyPlot.clf()
    # Print inhibitory assembly firing rates through time
    # gaus_kernel = ImageFiltering.Kernel.gaussian((fRate_window,))
    plot_data = zeros(Npop, T)
    for pp = 1:Npop
        plot_data[pp,:] = imfilter(iAssemfRates[pp,:], gaus_kernel)
    end
    labels = String[]
    xlim(0, T)
    ylim(0, maximum(plot_data))
    for pp = 1:Npop
        if mod(pp, seq_len) == 1 && pp > 1
            legend(labels)
            xlabel("Simulation time (ms)")
            ylabel("I₂-Assembly Firing Rate (Hz)")
            savefig(string("I-AssemfRate_", pp-seq_len, "-", pp-1, ".png"))
            PyPlot.clf()
            xlim(0, T)
            ylim(0, maximum(plot_data))
            labels = String[]
        end
        plot(x, plot_data[pp,:])
        append!(labels, [string("I-", pp)])
    end
    legend(labels)
    xlabel("Simulation time (ms)")
    ylabel("Assembly firing rate (Hz)")
    savefig(string("I-AssemfRate_", Npop-seq_len, "-", Npop, ".png"))
    PyPlot.clf()
end

function plotSingleTotalInput(Ncells, times, T, gaus_window=2)
    # ___________ Calculate total E & I input to a random neuron ___________
    # However, this does not include external excitatory input
    indxE = rand(1:Ne)                      # Get a random excitatory neuron
    indxi1 = rand((Ne+1):(Ncells-Ni2))      # Get a random neuron from 1st i-population
    indxi2 = rand((Ncells-Ni2+1):Ncells)    # Get a random neuron from 2nd i-population
    totalEinput = zeros(3, T)
    totalI1input = zeros(3, T)
    totalI2input = zeros(3, T)
    ii = 1
    for cc = 1:Ncells
        for tt in times[cc, :]
            (tt == 0) && continue
            (tt > 10) && (ii = Int(trunc(tt)))
            if cc <= Ne
                totalEinput[1, ii] += weights[cc, indxE]
                totalEinput[2, ii] += weights[cc, indxi1]
                totalEinput[3, ii] += weights[cc, indxi2]
            elseif cc <= (Ncells-Ni2)
                totalI1input[1, ii] += weights[cc, indxE]
                totalI1input[2, ii] += weights[cc, indxi1]
                totalI1input[3, ii] += weights[cc, indxi2]
            else
                totalI2input[1, ii] += weights[cc, indxE]
                totalI2input[2, ii] += weights[cc, indxi1]
                totalI2input[3, ii] += weights[cc, indxi2]
            end
        end
    end
    gaus_kernel = ImageFiltering.Kernel.gaussian((gaus_window,))
    y = 1:T
    # Plot for randomly chosen excitatory neuron
    plot(y, imfilter(totalEinput[1, :], gaus_kernel), c="blue")
    plot(y, imfilter(-totalI1input[1, :], gaus_kernel), c="magenta")
    plot(y, imfilter(-totalI2input[1, :], gaus_kernel), c="purple")
    plot(y, imfilter((-totalI2input[1, :]-totalI1input[1, :]), gaus_kernel), c="red")
    legend(["Excitatory", "1st i-population", "2nd i-population"])
    xlabel("Simulation time (ms)")
    ylabel(string("Input for excitatory neuron# ", indxE))
    savefig(string("totalInput_E.png"))
    PyPlot.clf()
    # Plot for randomly chosen 1st i-population neuron
    plot(y, imfilter(totalEinput[2, :], gaus_kernel), c="blue")
    plot(y, imfilter(-totalI1input[2, :], gaus_kernel), c="magenta")
    plot(y, imfilter(-totalI2input[2, :], gaus_kernel), c="purple")
    plot(y, imfilter((-totalI2input[2, :]-totalI1input[2, :]), gaus_kernel), c="red")
    legend(["Excitatory", "1st i-population", "2nd i-population"], fontsize="x-small")
    xlabel("Simulation time (ms)")
    ylabel(string("Input for 1st i-population neuron# ", indxi1))
    savefig(string("totalInput_i1.png"))
    PyPlot.clf()
    # Plot for randomly chosen 2nd i-population neuron
    plot(y, imfilter(totalEinput[3, :], gaus_kernel), c="blue")
    plot(y, imfilter(-totalI1input[3, :], gaus_kernel), c="magenta")
    plot(y, imfilter(-totalI2input[3, :], gaus_kernel), c="purple")
    plot(y, imfilter((-totalI2input[3, :]-totalI1input[3, :]), gaus_kernel), c="red")
    legend(["Excitatory", "1st i-population", "2nd i-population"], fontsize="small")
    xlabel("Simulation time (ms)")
    ylabel(string("Input for 2nd i-population neuron# ", indxi2))
    savefig(string("totalInput_i2.png"))
    PyPlot.clf()
end

# NOTE: NEW ADDITIONS TO THE CODE
Θ(x::Float64) = x > 0. ? x : 0.
function alpha_fun(t; t0, tau)
    (abs(t - t0)/ tau > 10) && (return 0.)
    return (t-t0) / tau * exp(1 - (t - t0) / tau) * Θ(1. * (t - t0))
end

# function convolveSpikes(spikeTimes::Vector{Float64}; interval::AbstractVector, tau=100)
#     # NOTE: HERE PROBABLY SPIKETIMES IS FOR A SINGLE NEURON (OR ACCUMULATED OF MANY IN ONE VECTOR)
#     # Taken from Alessio
#     rate = zeros(length(interval))
#     for t0 in spikeTimes
#         x = alpha_fun.(interval, t0=t0, tau=tau)
#         x[isnan.(x)] .= 0.
#         rate[:] .+= x
#     end
#     return rate
# end

# function convolveSpikes(spikeTimes; interval::AbstractVector, tau=100)
#     Npop = size(spikeTimes)[1]
#     rate = zeros(Npop, length(interval))
#     for ipop = 1:Npop
#         for t0 in filter(i->(i>0), spikeTimes[ipop, :])
#             x = alpha_fun.(interval, t0=t0, tau=tau)
#             x[isnan.(x)] .= 0.
#             rate[ipop, :] .+= x
#         end
#     end
#     return rate
# end

function getAssemblySpikes1(rates, popmembers, ipopmembers; Npop=20)
    T = size(rates)[2]
    eAssemfRates = zeros(Npop, T)
    iAssemfRates = zeros(Npop, T)
    for ipop = 1:Npop
        members = filter(i->(i>0), popmembers[ipop, :])
        eAssemfRates[ipop, :] = sum(rates[members, :], dims=1) ./ length(members)
        imembers = filter(i->(i>0), ipopmembers[ipop, :])
        iAssemfRates[ipop, :] = sum(rates[imembers, :], dims=1) ./ length(imembers)
    end
    restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
    inhfRates = sum(rates[4000:4500, :], dims=1) ./ 500
    restfRates = sum(rates[restcells, :], dims=1) ./ length(restcells)
    return eAssemfRates, iAssemfRates, restfRates, inhfRates
end

function getPopulationfRates(rates, Npop)
    T = size(rates)[2]
    emembers = sort(unique(popmembers))[2:end]
    restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), emembers)
    restfRates = sum(rates[restcells, :], dims=1) ./ length(restcells)
    excfRates = sum(rates[emembers, :], dims=1) ./ length(emembers)
    inhfRates = sum(rates[4001:4500, :], dims=1) ./ 500
    inh2fRates = sum(rates[4501:5000, :], dims=1) ./ 500
    return excfRates, restfRates, inhfRates, inh2fRates
end

function findAssemblyActivations(assemfRates; lowlim=0)
    # At each timestep return the most active assembly and if its firing rate is bigger than lowlim
    T = size(assemfRates)[2]
    succAssemAct = zeros(T)
    for tt = 1:T
        val, ind = findmax(assemfRates[:, tt])
        (val > lowlim) ? succAssemAct[tt] = ind : succAssemAct[tt] = 0.
    end
    return succAssemAct
end

function findSeqActivations(succAssemAct; seq_len=4, Npop=20)
    Nseq = round(Int, Npop/seq_len)
    succSeqAct = zeros(Nseq, length(succAssemAct))
    counter = zeros(Nseq, 2) # Contains
    for (tt, act) in enumerate(succAssemAct)
        (act == 0) && continue
        seq_ind = round(Int, div(act, seq_len) + 1) # Sequence number
        ass_ind = round(Int, mod(act, seq_len))     # Assembly number withing sequence
        if ass_ind == 0
            ass_ind = seq_len
            # (seq_ind > Nseq) ? seq_ind = Nseq : seq_ind -= 1
            seq_ind -= 1
        end
        if ass_ind == 1
            if counter[seq_ind, 2] != 1
                counter[seq_ind, 2] = 1
                counter[seq_ind, 1] = tt
                counter[1:Nseq .!= seq_ind, :] .= 0
            end
        elseif ass_ind == 2
            if counter[seq_ind, 2] != 2
                (counter[seq_ind, 2] == 1) ? counter[seq_ind, 2] = 2 : counter[seq_ind, :] .= 0
                counter[1:Nseq .!= seq_ind, :] .= 0
            end
        elseif ass_ind == 3
            if counter[seq_ind, 2] != 3
                (counter[seq_ind, 2] == 2) ? counter[seq_ind, 2] = 3 : counter[seq_ind, :] .= 0
                counter[1:Nseq .!= seq_ind, :] .= 0
            end
        else
            if counter[seq_ind, 2] != 4
                min_tt = round(Int, counter[seq_ind, 1])
                (counter[seq_ind, 2] == 3) ? (succSeqAct[seq_ind, min_tt:tt] .= 1) : counter[seq_ind, :] .= 0
                counter[1:Nseq .!= seq_ind, :] .= 0
            end
        end
    end
    return succSeqAct
end

function numSeqActivations(succSeqAct)
    # Returns the number of successful activations of whole sequences
    Nseq, T = size(succSeqAct)
    flags = zeros(Nseq)
    counts = zeros(Nseq)
    for tt = 1:T
        for iseq = 1:Nseq
            if flags[iseq] == 0 && succSeqAct[iseq, tt] == 1
                flags[iseq] = 1
            elseif flags[iseq] == 1 && succSeqAct[iseq, tt] == 0
                flags[iseq] = 0
                counts[iseq] += 1
            end
        end
    end
    return counts
end

function seqCorrelationCoef(succAssemAct; seq_len=4, Npop=20, interval=0)
    # Measures the percentage of feedforward dynamics activations
    # If activation of assembly i is followed by i+1
    if interval != 0
        succAssemAct = succAssemAct[interval]
    end
    Nseq = round(Int, Npop/seq_len)
    NassAct = zeros(Nseq)
    NseqAct = zeros(Nseq)
    prevAct = 0

    for act in succAssemAct
        if prevAct == act || act == 0
            continue
        else
            iSeq = ceil(Int, act/seq_len)
            if act in 1:seq_len:(Npop-seq_len+1)
                prevAct = act
            else
                if prevAct == (act-1)
                    NassAct[iSeq] += 1
                    NseqAct[iSeq] += 1
                    prevAct = act
                else
                    NassAct[iSeq] += 1
                    prevAct = act
                end
            end
        end
    end

    return (NseqAct, NassAct), sum(NseqAct)/sum(NassAct)
end

function seqEICorrelationCoef(succeAssemAct, succiAssemAct; seq_len=4, Npop=20, interval=0)
    # Measures the percentage of feedforward dynamics activations
    # If activation of assembly i is followed by i+1
    if interval != 0
        succeAssemAct = succeAssemAct[interval]
        succiAssemAct = succiAssemAct[interval]
    end
    Nseq = round(Int, Npop/seq_len)
    NasseAct = zeros(Nseq)
    NseqAct = zeros(Nseq)
    prevActe = 0
    for (tt, acte) in enumerate(succeAssemAct)
        if prevActe == acte || acte == 0 || tt == 1
            continue
        else
            iSeq = ceil(Int, acte/seq_len)
            acti = succiAssemAct[tt-1]
            if acte in 1:seq_len:(Npop-seq_len+1)
                prevActe = acte
            else
                if prevActe == (acte-1) && acti == (acte-1) # When feedforward dynamics succeed, check most active i Assembly
                    NasseAct[iSeq] += 1
                    NseqAct[iSeq] += 1
                    prevActe = acte
                elseif prevActe == (acte-1) # When feedforward dynamics succeed but the most active i Assembly is not the previous in order
                    NasseAct[iSeq] += 1
                    prevActe = acte
                else
                    prevActe = acte
                end
            end
        end
    end
    return (NseqAct, NasseAct), sum(NseqAct)/sum(NasseAct)
end


function calculateQ(succAssemAct; seq_len=4, Npop=20, interval=0)
    # Measures the network performance (parameter Q)
    # p_seq * (1 - variance(normalized(NassemAct)))
    if interval != 0
        succAssemAct = succAssemAct[interval]
    end

    current = 0
    expectation = 0
    NassemAct = zeros(Npop)
    NsuccAct = zeros(Npop)

    for act in succAssemAct # Calculate the number of successful sequential dynamics
        if current == act || act == 0
            continue
        else    # First moment of assembly activation
            NassemAct[round(Int, act)] += 1
            if act in 1:seq_len:(Npop-seq_len+1)
                expectation = act + 1
                current = act
            else
                if expectation == 0
                    current = act
                else
                    if act == expectation
                        NsuccAct[round(Int, act)] += 1
                        current = act
                    else
                        current = act
                    end
                end
                (act in seq_len:seq_len:Npop) ? expectation = 0 : expectation = act + 1
            end
        end
    end
    variance = var(NassemAct ./ maximum(NassemAct))
    deleteat!(NassemAct, 1:seq_len:(Npop-seq_len+1))
    deleteat!(NsuccAct, 1:seq_len:(Npop-seq_len+1))
    return (sum(NsuccAct)/sum(NassemAct)) * (1 - variance), variance
end

function plot_seqCorrelationCoef(succAssemAct, popmembers, ipopmembers; seq_len=4, Npop=20, interval=0)
    # Measures the percentage of feedforward dynamics activations
    # If activation of assembly i is followed by i+1
    if interval != 0
        succAssemAct = succAssemAct[interval]
    end
    Nseq = round(Int, Npop/seq_len)
    NsuccAct = zeros(Npop)
    NfailAct = zeros(Npop)
    current = 0
    expectation = 0
    NassemAct = zeros(Npop)

    for act in succAssemAct # Calculate the number of successful sequential dynamics
        if current == act || act == 0
            continue
        else    # First moment of assembly activation
            NassemAct[round(Int, act)] += 1
            if act in 1:seq_len:(Npop-seq_len+1)
                expectation = act + 1
                current = act
            else
                if expectation == 0
                    current = act
                else
                    if act == expectation
                        NsuccAct[round(Int, act)] += 1
                        current = act
                    else
                        NfailAct[round(Int, act)] += 1
                        current = act
                    end
                end
                (act in seq_len:seq_len:Npop) ? expectation = 0 : expectation = act + 1
            end
        end
    end
    # Nact = copy(NsuccAct)
    # Nact = copy(NfailAct)
    # Nact = copy(NfailAct-NsuccAct)
    # deleteat!(Nact, 1:seq_len:(Npop-seq_len+1))
    deleteat!(NfailAct, 1:seq_len:(Npop-seq_len+1))
    deleteat!(NsuccAct, 1:seq_len:(Npop-seq_len+1))
    # normalize!(Nact)

    dataEI = zeros(Npop)  # Contains the average strength from E-to-I2 assemblies
    for ipop = 1:Npop
        emembers = filter(i->(i>0), popmembers[ipop,:])
        imembers = ipopmembers[ipop, :]
        dataEI[ipop] = sum(new_weights[emembers, imembers]) / count(i->(i>0), new_weights[emembers, imembers])
    end
    # Make one plot here
    deleteat!(dataEI, seq_len:seq_len:Npop)
    # normalize!(dataEI)


    # dataIE = zeros(Npop)  # Contains the average strength from I2(i)-to-E(i+1)
    # for ipop = 2:Npop
    #     imembers = ipopmembers[ipop-1, :]
    #     emembers = filter(i->(i>0), popmembers[ipop, :])
    #     dataIE[ipop] = sum(new_weights[imembers, emembers]) / count(i->(i>0), new_weights[imembers, emembers])
    # end
    # # Make the other plot here
    # deleteat!(dataIE, 1:seq_len:(Npop-seq_len+1))
    # normalize!(dataIE)

    dataIE = zeros(Npop)  # Contains the average strength from I2(i)-to-E(i)
    for ipop = 1:Npop
        imembers = ipopmembers[ipop, :]
        emembers = filter(i->(i>0), popmembers[ipop, :])
        dataIE[ipop] = sum(new_weights[imembers, emembers]) / count(i->(i>0), new_weights[imembers, emembers])
        # dataIE[ipop] = sum(new_weights[4001:4500, emembers]) / count(i->(i>0), new_weights[4001:4500, emembers])
        # dataIE[ipop] = sum(new_weights[4001:5000, emembers]) / count(i->(i>0), new_weights[4001:5000, emembers])
    end
    # Make the other plot here
    deleteat!(dataIE, 1:seq_len:(Npop-seq_len+1))
    # normalize!(dataIE)

    summed = dataEI .+ dataIE
    # summed = normalize(dataIE) .+ -normalize(dataEI)
    # ______________________________

    # correlationEI = cor(dataEI, Nact)
    # correlationEI = corspearman(dataEI, Nact)
    # correlationEI = corkendall(dataEI, Nact)

    # correlationIE = cor(dataIE, Nact)
    # correlationIE = corspearman(dataIE, Nact)
    # correlationIE = corkendall(dataIE, Nact)

    # summedCor = cor(summed, Nact)
    # summedCor = corspearman(summed, Nact)
    # summedCor = corkendall(summed, Nact)
    # return correlationEI, correlationIE, summedCor
    return dataEI, dataIE, summed, NsuccAct, NfailAct, NassemAct
end

function plot_compare_eiSTDP(popmembers, ipopmembers, Npop, popmembers2, ipopmembers2, popmembers3, ipopmembers3)
    # THIS IS FOR I2-to-E weights
    avg_weights = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers[pp, :])
        imembers = ipopmembers[pp, :]
        avg_weights[pp] = sum(new_weights[imembers, emembers]) / count(i->i>0,new_weights[imembers, emembers])
    end
    avg_weights2 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers2[pp, :])
        imembers = ipopmembers2[pp, :]
        avg_weights2[pp] = sum(w8s[imembers, emembers]) / count(i->i>0,w8s[imembers, emembers])
    end
    avg_weights3 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers3[pp, :])
        imembers = ipopmembers3[pp, :]
        avg_weights3[pp] = sum(w8s2[imembers, emembers]) / count(i->i>0,w8s2[imembers, emembers])
    end
    avg_weights4 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers4[pp, :])
        imembers = ipopmembers4[pp, :]
        avg_weights4[pp] = sum(w8s3[imembers, emembers]) / count(i->i>0,w8s3[imembers, emembers])
    end
    # x = 1:Npop
    figure(figsize=(12, 4))
    # xlim(0, Npop+1)
    # ylim(0, 360)
    plt = matplotlib.pyplot.boxplot([avg_weights, avg_weights2, avg_weights4, avg_weights3], labels=["with\neiSTDP", "without\neiSTDP", "with\neiSTDP\n(no regulation)", "without\neiSTDP\n(no regulation)"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2, "color"=>"tab:orange"))
    [box.set_facecolor("tab:red") for box in plt["boxes"]]
    # boxplot([avg_weights, avg_weights3, avg_weights2],labels=["with\neiSTDP", "without\neiSTDP\n(no regulation)", "without\neiSTDP"])
    ax = gca()
    ylabel("Average I₂→E Assembly\nCoupling Strength (pF)", fontsize="xx-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
    yticks(fontsize="x-large")
    savefig("boxI2-E_strength.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_compare_tau_IE()
    # THIS IS FOR I2-to-E weights
    avg_weights20 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers20[pp, :])
        imembers = ipopmembers20[pp, :]
        avg_weights20[pp] = sum(w8s20[imembers, emembers]) / count(i->i>0, w8s20[imembers, emembers])
    end
    avg_weights50 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers50[pp, :])
        imembers = ipopmembers50[pp, :]
        avg_weights50[pp] = sum(w8s50[imembers, emembers]) / count(i->i>0, w8s50[imembers, emembers])
    end
    avg_weights100 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers100[pp, :])
        imembers = ipopmembers100[pp, :]
        avg_weights100[pp] = sum(w8s100[imembers, emembers]) / count(i->i>0, w8s100[imembers, emembers])
    end
    avg_weights200 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers200[pp, :])
        imembers = ipopmembers200[pp, :]
        avg_weights200[pp] = sum(w8s200[imembers, emembers]) / count(i->i>0, w8s200[imembers, emembers])
    end
    avg_weights500 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers500[pp, :])
        imembers = ipopmembers500[pp, :]
        avg_weights500[pp] = sum(w8s500[imembers, emembers]) / count(i->i>0, w8s500[imembers, emembers])
    end
    avg_weights20 = mean(avg_weights20, dims=2)
    avg_weights50 = mean(avg_weights50, dims=2)
    avg_weights100 = mean(avg_weights100, dims=2)
    avg_weights200 = mean(avg_weights200, dims=2)
    avg_weights500 = mean(avg_weights500, dims=2)
    figure(figsize=(8, 4))
    plt = boxplot([avg_weights20, avg_weights50, avg_weights100, avg_weights200, avg_weights500],labels=["20", "50", "100", "200", "500"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2.5, "color"=>"tab:cyan"))
    [box.set_facecolor("tab:red") for box in plt["boxes"]]
    ax = gca()
    xlabel("iSTDP₂ time constant (ms)", fontsize="x-large")
    ylabel("Average I₂→E Assembly\nCoupling Strength (pF)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    yticks(fontsize="large")
    savefig("boxTauI2-E.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_compare_tau_EI()
    # THIS IS FOR E-to-I2 weights
    avg_weights20 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers20[pp, :])
        imembers = ipopmembers20[pp, :]
        avg_weights20[pp] = sum(w8s20[emembers, imembers]) / count(i->i>0, w8s20[emembers, imembers])
    end
    avg_weights50 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers50[pp, :])
        imembers = ipopmembers50[pp, :]
        avg_weights50[pp] = sum(w8s50[emembers, imembers]) / count(i->i>0, w8s50[emembers, imembers])
    end
    avg_weights100 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers100[pp, :])
        imembers = ipopmembers100[pp, :]
        avg_weights100[pp] = sum(w8s100[emembers, imembers]) / count(i->i>0, w8s100[emembers, imembers])
    end
    avg_weights200 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers200[pp, :])
        imembers = ipopmembers200[pp, :]
        avg_weights200[pp] = sum(w8s200[emembers, imembers]) / count(i->i>0, w8s200[emembers, imembers])
    end
    avg_weights500 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers500[pp, :])
        imembers = ipopmembers500[pp, :]
        avg_weights500[pp] = sum(w8s500[emembers, imembers]) / count(i->i>0, w8s500[emembers, imembers])
    end
    avg_weights20 = mean(avg_weights20, dims=2)
    avg_weights50 = mean(avg_weights50, dims=2)
    avg_weights100 = mean(avg_weights100, dims=2)
    avg_weights200 = mean(avg_weights200, dims=2)
    avg_weights500 = mean(avg_weights500, dims=2)
    figure(figsize=(8, 4))
    plt = boxplot([avg_weights20, avg_weights50, avg_weights100, avg_weights200, avg_weights500],labels=["20", "50", "100", "200", "500"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2.5, "color"=>"tab:orange"))
    [box.set_facecolor("tab:blue") for box in plt["boxes"]]
    ax = gca()
    xlabel("iSTDP₂ time constant (ms)", fontsize="x-large")
    ylabel("Average E→I₂ Assembly\nCoupling Strength (pF)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    yticks(fontsize="large")
    savefig("boxTauE-I2.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_compare_ilamda_IE()
    # THIS IS FOR I2-to-E weights
    avg_weights01 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers01[pp, :])
        imembers = ipopmembers01[pp, :]
        avg_weights01[pp] = sum(w8s01[imembers, emembers]) / count(i->i>0, w8s01[imembers, emembers])
    end
    avg_weights05 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers05[pp, :])
        imembers = ipopmembers05[pp, :]
        avg_weights05[pp] = sum(w8s05[imembers, emembers]) / count(i->i>0, w8s05[imembers, emembers])
    end
    avg_weights1 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers1[pp, :])
        imembers = ipopmembers1[pp, :]
        avg_weights1[pp] = sum(w8s1[imembers, emembers]) / count(i->i>0, w8s1[imembers, emembers])
    end
    avg_weights10 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers10[pp, :])
        imembers = ipopmembers10[pp, :]
        avg_weights10[pp] = sum(w8s10[imembers, emembers]) / count(i->i>0, w8s10[imembers, emembers])
    end
    avg_weights100 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers100[pp, :])
        imembers = ipopmembers100[pp, :]
        avg_weights100[pp] = sum(w8s100[imembers, emembers]) / count(i->i>0, w8s100[imembers, emembers])
    end
    avg_weights500 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers500[pp, :])
        imembers = ipopmembers500[pp, :]
        avg_weights500[pp] = sum(w8s500[imembers, emembers]) / count(i->i>0, w8s500[imembers, emembers])
    end
    avg_weights01 = mean(avg_weights01, dims=2)
    avg_weights05 = mean(avg_weights05, dims=2)
    avg_weights1 = mean(avg_weights1, dims=2)
    avg_weights10 = mean(avg_weights10, dims=2)
    avg_weights100 = mean(avg_weights100, dims=2)
    avg_weights500 = mean(avg_weights500, dims=2)
    figure(figsize=(8, 4))
    plt = boxplot([avg_weights01, avg_weights05, avg_weights1, avg_weights10, avg_weights100, avg_weights500],labels=["01", "05", "1", "10", "100", "500"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2.5, "color"=>"tab:cyan"))
    [box.set_facecolor("tab:red") for box in plt["boxes"]]
    ax = gca()
    xlabel("iSTDP₂ learning rate", fontsize="x-large")
    ylabel("Average I₂→E Assembly\nCoupling Strength (pF)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    yticks(fontsize="large")
    savefig("boxTauI2-E.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_compare_ilamda_EI()
    # THIS IS FOR E-to-I2 weights
    avg_weights01 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers01[pp, :])
        imembers = ipopmembers01[pp, :]
        avg_weights01[pp] = sum(w8s01[emembers, imembers]) / count(i->i>0, w8s01[emembers, imembers])
    end
    avg_weights05 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers05[pp, :])
        imembers = ipopmembers05[pp, :]
        avg_weights05[pp] = sum(w8s05[emembers, imembers]) / count(i->i>0, w8s05[emembers, imembers])
    end
    avg_weights1 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers1[pp, :])
        imembers = ipopmembers1[pp, :]
        avg_weights1[pp] = sum(w8s1[emembers, imembers]) / count(i->i>0, w8s1[emembers, imembers])
    end
    avg_weights10 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers10[pp, :])
        imembers = ipopmembers10[pp, :]
        avg_weights10[pp] = sum(w8s10[emembers, imembers]) / count(i->i>0, w8s10[emembers, imembers])
    end
    avg_weights100 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers100[pp, :])
        imembers = ipopmembers100[pp, :]
        avg_weights100[pp] = sum(w8s100[emembers, imembers]) / count(i->i>0, w8s100[emembers, imembers])
    end
    avg_weights500 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers500[pp, :])
        imembers = ipopmembers500[pp, :]
        avg_weights500[pp] = sum(w8s500[emembers, imembers]) / count(i->i>0, w8s500[emembers, imembers])
    end
    avg_weights01 = mean(avg_weights01, dims=2)
    avg_weights05 = mean(avg_weights05, dims=2)
    avg_weights1 = mean(avg_weights1, dims=2)
    avg_weights10 = mean(avg_weights10, dims=2)
    avg_weights100 = mean(avg_weights100, dims=2)
    avg_weights500 = mean(avg_weights500, dims=2)
    figure(figsize=(8, 4))
    plt = boxplot([avg_weights01, avg_weights05, avg_weights1, avg_weights10, avg_weights100, avg_weights500],labels=["01", "05", "1", "10", "100", "500"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2.5, "color"=>"tab:orange"))
    [box.set_facecolor("tab:blue") for box in plt["boxes"]]
    ax = gca()
    xlabel("iSTDP₂ learning rate", fontsize="x-large")
    ylabel("Average E→I₂ Assembly\nCoupling Strength (pF)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    yticks(fontsize="large")
    savefig("boxTauE-I2.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_compare_ii_IE()
    # THIS IS FOR I2-to-E weights
    avg_weights16 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers16[pp, :])
        imembers = ipopmembers16[pp, :]
        avg_weights16[pp] = sum(weights16[imembers, emembers]) / count(i->i>0, weights16[imembers, emembers])
    end
    avg_weights20 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers20[pp, :])
        imembers = ipopmembers20[pp, :]
        avg_weights20[pp] = sum(weights20[imembers, emembers]) / count(i->i>0, weights20[imembers, emembers])
    end
    avg_weights28 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers28[pp, :])
        imembers = ipopmembers28[pp, :]
        avg_weights28[pp] = sum(weights28[imembers, emembers]) / count(i->i>0, weights28[imembers, emembers])
    end
    avg_weightsOrig = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers[pp, :])
        imembers = ipopmembers[pp, :]
        avg_weightsOrig[pp] = sum(weightsOrig[imembers, emembers]) / count(i->i>0, weightsOrig[imembers, emembers])
    end
    avg_weights32 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers32[pp, :])
        imembers = ipopmembers32[pp, :]
        avg_weights32[pp] = sum(weights32[imembers, emembers]) / count(i->i>0, weights32[imembers, emembers])
    end
    avg_weights40 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers40[pp, :])
        imembers = ipopmembers40[pp, :]
        avg_weights40[pp] = sum(weights40[imembers, emembers]) / count(i->i>0, weights40[imembers, emembers])
    end
    avg_weights50 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers50[pp, :])
        imembers = ipopmembers50[pp, :]
        avg_weights50[pp] = sum(weights50[imembers, emembers]) / count(i->i>0, weights50[imembers, emembers])
    end
    avg_weights16 = mean(avg_weights16, dims=2)
    avg_weights20 = mean(avg_weights20, dims=2)
    avg_weights28 = mean(avg_weights28, dims=2)
    avg_weightsOrig = mean(avg_weightsOrig, dims=2)
    avg_weights32 = mean(avg_weights32, dims=2)
    avg_weights40 = mean(avg_weights40, dims=2)
    avg_weights50 = mean(avg_weights50, dims=2)
    figure(figsize=(8, 4))
    plt = boxplot([avg_weights16, avg_weights20, avg_weights28, avg_weightsOrig, avg_weights32, avg_weights40, avg_weights50],labels=["16", "20", "28", "30", "32", "40", "50"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2.5, "color"=>"tab:cyan"))
    [box.set_facecolor("tab:red") for box in plt["boxes"]]
    ax = gca()
    xlabel("I₁→I₂ Synaptic Strength (pF)", fontsize="x-large")
    ylabel("Average I₂→E Assembly\nCoupling Strength (pF)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    yticks(fontsize="large")
    savefig("boxiiI2-E.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_compare_ii_EI()
    # THIS IS FOR I2-to-E weights
    avg_weights16 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers16[pp, :])
        imembers = ipopmembers16[pp, :]
        avg_weights16[pp] = sum(weights16[emembers, imembers]) / count(i->i>0, weights16[emembers, imembers])
    end
    avg_weights20 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers20[pp, :])
        imembers = ipopmembers20[pp, :]
        avg_weights20[pp] = sum(weights20[emembers, imembers]) / count(i->i>0, weights20[emembers, imembers])
    end
    avg_weights28 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers28[pp, :])
        imembers = ipopmembers28[pp, :]
        avg_weights28[pp] = sum(weights28[emembers, imembers]) / count(i->i>0, weights28[emembers, imembers])
    end
    avg_weightsOrig = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembersOrig[pp, :])
        imembers = ipopmembersOrig[pp, :]
        avg_weightsOrig[pp] = sum(weightsOrig[emembers, imembers]) / count(i->i>0, weightsOrig[emembers, imembers])
    end
    avg_weights32 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers32[pp, :])
        imembers = ipopmembers32[pp, :]
        avg_weights32[pp] = sum(weights32[emembers, imembers]) / count(i->i>0, weights32[emembers, imembers])
    end
    avg_weights40 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers40[pp, :])
        imembers = ipopmembers40[pp, :]
        avg_weights40[pp] = sum(weights40[emembers, imembers]) / count(i->i>0, weights40[emembers, imembers])
    end
    avg_weights50 = zeros(Npop)
    for pp = 1:Npop
        emembers = filter(i->(i>0), popmembers50[pp, :])
        imembers = ipopmembers50[pp, :]
        avg_weights50[pp] = sum(weights50[emembers, imembers]) / count(i->i>0, weights50[emembers, imembers])
    end
    avg_weights16 = mean(avg_weights16, dims=2)
    avg_weights20 = mean(avg_weights20, dims=2)
    avg_weights28 = mean(avg_weights28, dims=2)
    avg_weightsOrig = mean(avg_weightsOrig, dims=2)
    avg_weights32 = mean(avg_weights32, dims=2)
    avg_weights40 = mean(avg_weights40, dims=2)
    avg_weights50 = mean(avg_weights50, dims=2)
    figure(figsize=(8, 4))
    plt = boxplot([avg_weights16, avg_weights20, avg_weights28, avg_weightsOrig, avg_weights32, avg_weights40, avg_weights50],labels=["16", "20", "28", "30", "32", "40", "50"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2.5, "color"=>"tab:orange"))
    [box.set_facecolor("tab:blue") for box in plt["boxes"]]
    ax = gca()
    xlabel("I₁→I₂ Synaptic Strength (pF)", fontsize="x-large")
    ylabel("Average E→I₂ Assembly\nCoupling Strength (pF)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    yticks(fontsize="large")
    savefig("boxiiE-I2.png", bbox_inches="tight")
    PyPlot.clf()
end

function measureInterference(weights, popmembers, ipopmembers; Npop=20)
    totalInputI1 = zeros(Npop)
    averageInputI1 = zeros(Npop)
    totalInputI2 = zeros(Npop)
    averageInputI2 = zeros(Npop)
    for ipop = 1:Npop
        emembers = filter(i->(i>0), popmembers[ipop, :])
        imembers = ipopmembers[ipop, :]
        totalInputI2[ipop] = sum(weights[imembers, emembers])
        averageInputI2[ipop] = totalInputI2[ipop] / count(i->(i>0), weights[imembers, emembers])
        totalInputI1[ipop] = sum(weights[(Ne+1):(Ncells-Ni2), emembers])
        averageInputI1[ipop] = totalInputI1[ipop] / count(i->(i>0), weights[(Ne+1):(Ncells-Ni2), emembers])
    end
    return mean(averageInputI2 ./ averageInputI1), averageInputI2 ./ averageInputI1#mean(totalInputI2 ./ totalInputI1)
end

function calculateEcoupling(weights, popmembers; Npop=20)
    populations = length(popmembers)
    avgE = zeros(populations, Npop)
    for pops = 1:populations
        for ipop = 1:Npop
            emembers = filter(i->(i>0), popmembers[pops][ipop, :])
            # avgE[pops, ipop] = mean(weights[pops][emembers, emembers])
            avgE[pops, ipop] = sum(weights[pops][emembers, emembers]) / count(i->(i>0), weights[pops][emembers, emembers])
        end
    end
    return avgE
end

function boxplot_IIinter()
    figure(figsize=(12, 4))
    # xlim(0, Npop+1)
    ylim(0, 19)
    plt = matplotlib.pyplot.boxplot(avgEcoupling', labels=["16", "20", "28", "30", "32", "40", "50"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2, "color"=>"tab:orange"))
    [box.set_facecolor("tab:blue") for box in plt["boxes"]]
    # [total_inter16, total_inter20, total_inter28, total_interOrig, total_inter32, total_inter40, total_inter50]
    ax = gca()
    xlabel("I₁→I₂ Synaptic Strength (pF)", fontsize="xx-large")
    ylabel("Average E Assembly\nCoupling Strength (pF)", fontsize="xx-large")
    ax.yaxis.set_tick_params(which="major", labelleft=false, labelright=true, left=false, right=true, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
    yticks(fontsize="x-large")
    plot(1:7, [avg_inter16, avg_inter20, avg_inter28, avg_interOrig, avg_inter32, avg_inter40, avg_inter50], label="Interference", color="tab:green", linewidth=2, marker="D")
    grid(true, axis="y")
    ax.text(0.6, 3, "Interference (Dᵢ)", color="tab:gray", fontsize="xx-large")
    savefig("boxII_interference.png", bbox_inches="tight")
    PyPlot.clf()
end

function plotCorrelation(data1, data2)
    fit = curve_fit(LinearFit, data1, data2)
    x = minimum(data1):((maximum(data1)-minimum(data1))/10):maximum(data1)
    y = fit.(x)
    figure(figsize=(8, 4))
    scatter(data1, data2, c="tab:red")
    plot(x, y, linewidth=3)
    ax = gca()
    # xlabel("Average I₂→E Assembly Coupling Strength", fontsize="x-large")
    # ylabel("Number of Failed Transitions", fontsize="x-large")
    xlabel("Total Inhibitory Input to E-assembly", fontsize="x-large")
    ylabel("E-assembly Activations", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    yticks(fontsize="large")
    props = Dict("boxstyle"=>"round", "facecolor"=>"tab:red", "alpha"=>0.5)
    # boxstyle='round', facecolor='wheat', alpha=0.5
    ax.text(0.05, 0.95, string("ρ: ", round(cor(data1, data2), digits=2)), transform=ax.transAxes, fontsize=14, verticalalignment="top", bbox=props)
    savefig("correlation.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_PopulationfRates()
    # Tau
    figure(figsize=(8, 4))
    # xlim(0.5,5.5)
    ylim(0, 3.5)
    exc = [mean(excfRates20), mean(excfRates50), mean(excfRates100), mean(excfRates200), mean(excfRates500)]
    inh1 = [mean(inhfRates20), mean(inhfRates50), mean(inhfRates100), mean(inhfRates200), mean(inhfRates500)]
    inh2 = [mean(inh2fRates20), mean(inh2fRates50), mean(inh2fRates100), mean(inh2fRates200), mean(inh2fRates500)]
    rest = [mean(restfRates20), mean(restfRates50), mean(restfRates100), mean(restfRates200), mean(restfRates500)]
    x = 1:5
    plot(x, exc, color="tab:blue", label="Excitatory (members)", linewidth=2.5, marker="D")
    plot(x, rest, color="tab:gray", label="Excitatory (non-members)", linewidth=2, marker="D", zorder=1)
    plot(x, inh1, color="brown", label="Inhibitory₁", linewidth=2, marker="D")
    plot(x, inh2, color="red", label="Inhibitory₂", linewidth=2.5, marker="D")
    ax = gca()
    xticks = matplotlib.ticker.FixedLocator(0:6)
	ax.xaxis.set_major_locator(xticks)
    xlabel("iSTDP₂ time constant (ms)", fontsize="x-large")
    ylabel("Mean Firing Rate (Hz)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    ax.set_xticklabels(["", "20", "50", "100", "200", "500", ""], minor=false, fontsize="large")
    yticks(fontsize="large")
    grid(true, axis="x")
    grid(true, axis="y", linestyle="--")
    lgd = legend(loc="upper left", ncol=2, bbox_to_anchor=(0,1), fontsize="large")
    savefig("boxfRates.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_PopulationfRates()
    # ilamda
    figure(figsize=(8, 4))
    # xlim(0.5,5.5)
    ylim(0, 3.5)
    exc = [mean(excfRates01), mean(excfRates05), mean(excfRates1), mean(excfRates10), mean(excfRates100), mean(excfRates500)]
    inh1 = [mean(inhfRates01), mean(inhfRates05), mean(inhfRates1), mean(inhfRates10), mean(inhfRates100), mean(inhfRates500)]
    inh2 = [mean(inh2fRates01), mean(inh2fRates05), mean(inh2fRates1), mean(inh2fRates10), mean(inh2fRates100), mean(inh2fRates500)]
    rest = [mean(restfRates01), mean(restfRates05), mean(restfRates1), mean(restfRates10), mean(restfRates100), mean(restfRates500)]
    x = 1:6
    plot(x, exc, color="tab:blue", label="Excitatory (members)", linewidth=2.5, marker="D")
    plot(x, rest, color="tab:gray", label="Excitatory (non-members)", linewidth=2, marker="D", zorder=1)
    plot(x, inh1, color="brown", label="Inhibitory₁", linewidth=2, marker="D")
    plot(x, inh2, color="red", label="Inhibitory₂", linewidth=2.5, marker="D")
    ax = gca()
    xticks = matplotlib.ticker.FixedLocator(0:7)
	ax.xaxis.set_major_locator(xticks)
    xlabel("iSTDP₂ learning rate", fontsize="x-large")
    ylabel("Mean Firing Rate (Hz)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    ax.set_xticklabels(["", "0.1", "0.5", "1", "10", "100", "500", ""], minor=false, fontsize="large")
    yticks(fontsize="large")
    grid(true, axis="x")
    grid(true, axis="y", linestyle="--")
    lgd = legend(loc="upper left", ncol=2, bbox_to_anchor=(0,1), fontsize="large")
    savefig("boxfRates.png", bbox_inches="tight")
    PyPlot.clf()
end

function plot_iPopulationfRates()
    # i1i2
    figure(figsize=(8, 4))
    # xlim(0.5,5.5)
    ylim(0, 3.2)
    exc = [y16[1], y20[1], y28[1], yOrig[1], y32[1], y40[1], y50[1]]
    inh1 = [y16[4], y20[4], y28[4], yOrig[4], y32[4], y40[4], y50[4]]
    inh2 = [y16[2], y20[2], y28[2], yOrig[2], y32[2], y40[2], y50[2]]
    rest = [y16[3], y20[3], y28[3], yOrig[3], y32[3], y40[3], y50[3]]
    x = 1:7
    plot(x, exc, color="tab:blue", label="Excitatory (members)", linewidth=2.5, marker="D")
    plot(x, rest, color="tab:gray", label="Excitatory (non-members)", linewidth=2, marker="D", zorder=1)
    plot(x, inh1, color="brown", label="Inhibitory₁", linewidth=2, marker="D")
    plot(x, inh2, color="red", label="Inhibitory₂", linewidth=2.5, marker="D")
    ax = gca()
    xticks = matplotlib.ticker.FixedLocator(0:7)
	ax.xaxis.set_major_locator(xticks)
    xlabel("I₁→I₂ Synaptic Strength (pF)", fontsize="x-large")
    ylabel("Mean Firing Rate (Hz)", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="large")
    ax.set_xticklabels(["", "16", "20", "28", "30", "32", "40", "50"], minor=false, fontsize="large")
    yticks(fontsize="large")
    grid(true, axis="x")
    grid(true, axis="y", linestyle="--")
    lgd = legend(loc="upper left", ncol=2, bbox_to_anchor=(0,1), fontsize="large")
    savefig("boxiifRates.png", bbox_inches="tight")
    PyPlot.clf()
end

function plotGroupfRates(excRates, inh1Rates, inh2Rates, restRates; mint=1000, maxt=6000)
    figure(figsize=(12, 3))
    xlim(0, maxt-mint)
    ylim(0, 3.1)
    ylabel("Firing Rate (Hz)", fontsize="large")
    xlabel("Simulation Time (ms)", fontsize="large")
    tight_layout()
    plot(excRates[mint:maxt], color="tab:blue", label="Excitatory (members)")
    plot(restRates[mint:maxt], color="tab:gray", label="Excitatory (non-members)")
    plot(inh1Rates[mint:maxt], color="brown", label="Inhibitory₁")
    plot(inh2Rates[mint:maxt], color="red", label="Inhibitory₂")
    ax = gca()
    lgd = legend(loc="upper left", ncol=4, bbox_to_anchor=(0,1), fontsize="large")
    xticks = matplotlib.ticker.FixedLocator(0:500:(maxt-mint))
	ax.xaxis.set_major_locator(xticks)
    # ax.set_xticklabels(0:500:(maxt-mint), minor=false)
    savefig(string("fireRate", T, ".png"))
    PyPlot.clf()
end

function plotGroupfRates2(excRates, inh1Rates, inh2Rates, restRates)
    figure(figsize=(10.5, 3))
    xlim(0, length(excfRates))
    ylim(0, 3.8)
    lineWidth = 1
    ylabel("Firing Rate (Hz)", fontsize="xx-large")
    # xlabel("Simulation Time (ms)", fontsize="large")
    tight_layout()
    plot(excRates, color="tab:blue", label="Excitatory (members)", lw=lineWidth)
    plot(restRates, color="tab:gray", label="Excitatory (non-members)", lw=lineWidth)
    plot(inh1Rates, color="brown", label="Inhibitory₁", lw=lineWidth)
    plot(inh2Rates, color="red", label="Inhibitory₂", lw=lineWidth)
    ax = gca()
    legend(loc="upper left", ncol=4, bbox_to_anchor=(0,1), fontsize="large")
    # xticks = matplotlib.ticker.FixedLocator(0:500:5000)
	# ax.xaxis.set_major_locator(xticks)
    ax.xaxis.set_tick_params(labeltop=false, labelbottom=false, top=false, bottom=true, labelsize="x-large")
    ax.yaxis.set_tick_params(labelright=false, labelleft=true, right=false, left=true, labelsize="x-large")
    savefig(string("fireRate", T, ".png"), bbox_inches="tight", dpi=300)
    PyPlot.clf()
end

function plotGroupfRates(excRates, inh1Rates, inh2Rates; mint=1000, maxt=6000)
    figure(figsize=(10.5, 2.3))
    xlim(0, maxt-mint)
    ylim(0, 5.2)
    lineWidth = 1.3  # .8
    ylabel("Average\nFiring Rate (Hz)", fontsize=28, labelpad=9)
    # xlabel("Simulation Time (ms)", fontsize=21)
    tight_layout()
    plot(excRates[mint:maxt], color="tab:blue", label="Excitatory", lw=lineWidth)
    # plot(restRates[mint:maxt], color="tab:gray", label="Excitatory (non-members)")
    plot(inh1Rates[mint:maxt], color="brown", label="Inhibitory I", lw=lineWidth)
    plot(inh2Rates[mint:maxt], color="red", label="Inhibitory II", lw=lineWidth)
    ax = gca()
    lgd = legend(loc="upper left", ncol=4, bbox_to_anchor=(0,1), fontsize="xx-large")
    # x_ticks = matplotlib.ticker.FixedLocator(0:1000:(maxt-mint))
	# ax.xaxis.set_major_locator(x_ticks)
    ax.xaxis.set_tick_params(labeltop=false, labelbottom=false, top=false, bottom=true, labelsize="xx-large")
    y_ticks = matplotlib.ticker.FixedLocator(0:2:4)
    ax.yaxis.set_major_locator(y_ticks)
    yticks(fontsize="xx-large")
    savefig(string("fireRate", T, ".png"), bbox_inches="tight", dpi=150)
    PyPlot.clf()
end

function plotGroupfRates(excRates, inhRates; mint=1000, maxt=6000)
    figure(figsize=(10.5, 2.3))
    xlim(0, maxt-mint)
    ylim(0, 4.3)
    lineWidth = 1.3 # .8
    ylabel("Average\nFiring Rate (Hz)", fontsize=28, labelpad=9)
    # xlabel("Simulation Time (ms)", fontsize="large")
    tight_layout()
    plot(excRates[mint:maxt], color="tab:blue", label="Excitatory", lw=lineWidth)
    plot(inhRates[mint:maxt], color="brown", label="Inhibitory", lw=lineWidth)    # brown/red
    ax = gca()
    lgd = legend(loc="upper left", ncol=4, bbox_to_anchor=(0,1), fontsize="xx-large")
    ax.xaxis.set_tick_params(labeltop=false, labelbottom=false, top=false, bottom=true, labelsize="xx-large")
    x_ticks = matplotlib.ticker.FixedLocator(0:1000:(maxt-mint))
	ax.xaxis.set_major_locator(x_ticks)
    ax.xaxis.set_tick_params(labelsize="xx-large")
    yticks(fontsize="xx-large")
    savefig(string("fireRate", T, ".png"), bbox_inches="tight", dpi=150)
    PyPlot.clf()
end

function plotWeightDist_tau()
    figure(figsize=(8, 6))
    xlim(28,210)
    ylim(0, 0.08)
    ylabel("Frequency", fontsize="xx-large")
    xlabel("Average I₂→E Assembly Coupling Strength (pF)", fontsize="xx-large")
    nbins = 60
    alpha = 0.3 # opacity
    ln_width = 3.5
    flt_window = 2.5
    hcolors = matplotlib.cm.Reds_r(1:40:200)
    pcolors = copy(hcolors)
    hcolors[:, 4] .= alpha
    # τ20
    n, bins, patches = hist(x20, nbins, density=true, fc=hcolors[1, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[1, :], label=string("τ=20,   Q=", round(q20, digits=2), ", (μ=", round(mean(x20), digits=1), ", σ=", round(std(x20), digits=1), ")"), linewidth=ln_width)
    # τ50
    n, bins, patches = hist(x50, nbins, density=true, fc=hcolors[2, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[2, :], label=string("τ=50,   Q=", round(q50, digits=2), ", (μ=", round(mean(x50), digits=1), ", σ=", round(std(x50), digits=1), ")"), linewidth=ln_width)
    # τ100
    n, bins, patches = hist(xOrig, nbins, density=true, fc=hcolors[3, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[3, :], label=string("*τ=100, Q=", round(qOrig, digits=2), ", (μ=", round(mean(xOrig), digits=1), ", σ=", round(std(xOrig), digits=1), ")"), linewidth=ln_width)
    # τ200
    n, bins, patches = hist(x200, nbins, density=true, fc=hcolors[4, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[4, :], label=string("τ=200, Q=", round(q200, digits=2), ", (μ=", round(mean(x200), digits=1), ", σ=", round(std(x200), digits=1), ")"), linewidth=ln_width)
    # τ500
    n, bins, patches = hist(x500, nbins, density=true, fc=hcolors[5, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[5, :], label=string("τ=500, Q=", round(q500_t, digits=2), ", (μ=", round(mean(x500), digits=1), ", σ=", round(std(x500), digits=1), ")"), linewidth=ln_width)

    ax = gca()
    lgd = legend(loc="upper right", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
    savefig(string("weight_dist_tau.png"), bbox_inches="tight")
    PyPlot.clf()
end

function plotWeightDist_ilamda()
    figure(figsize=(8, 6))
    xlim(28,450)
    ylim(0, 0.09)
    ylabel("Frequency", fontsize="xx-large")
    xlabel("Average I₂→E Assembly Coupling Strength (pF)", fontsize="xx-large")
    nbins = 60
    alpha = 0.3 # opacity
    ln_width = 3.5
    flt_window = 2.5
    hcolors = matplotlib.cm.Reds_r(1:40:200)
    pcolors = copy(hcolors)
    hcolors[:, 4] .= alpha
    # τ0.1
    n, bins, patches = hist(x01, nbins, density=true, fc=hcolors[1, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[1, :], label=string("λ=0.1,  Q=", round(q01, digits=2), ", (μ=", round(mean(x01), digits=1), ", σ=", round(std(x01), digits=1), ")"), linewidth=ln_width)
    # τ0.5
    n, bins, patches = hist(x05, nbins, density=true, fc=hcolors[2, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[2, :], label=string("λ=0.5,  Q=", round(q05, digits=2), ", (μ=", round(mean(x05), digits=1), ", σ=", round(std(x05), digits=1), ")"), linewidth=ln_width)
    # τ10
    n, bins, patches = hist(xOrig, nbins, density=true, fc=hcolors[3, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[3, :], label=string("*λ=1, Q=", round(qOrig, digits=2), ", (μ=", round(mean(xOrig), digits=1), ", σ=", round(std(xOrig), digits=1), ")"), linewidth=ln_width)
    # τ100
    n, bins, patches = hist(x10, nbins, density=true, fc=hcolors[3, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[3, :], label=string("λ=10, Q=", round(q10, digits=2), ", (μ=", round(mean(x10), digits=1), ", σ=", round(std(x10), digits=1), ")"), linewidth=ln_width)
    # τ200
    n, bins, patches = hist(x100, nbins, density=true, fc=hcolors[4, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[4, :], label=string("λ=100, Q=", round(q100, digits=2), ", (μ=", round(mean(x100), digits=1), ", σ=", round(std(x100), digits=1), ")"), linewidth=ln_width)
    # τ500
    n, bins, patches = hist(x500_i, nbins, density=true, fc=hcolors[5, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[5, :], label=string("λ=500, Q=", round(q500_i, digits=2), ", (μ=", round(mean(x500_i), digits=1), ", σ=", round(std(x500_i), digits=1), ")"), linewidth=ln_width)

    ax = gca()
    lgd = legend(loc="upper right", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
    savefig(string("weight_dist_ilamda.png"), bbox_inches="tight")
    PyPlot.clf()
end

function plotWeightDist_ii()
    figure(figsize=(8, 6))
    xlim(28, 150)
    ylim(0, 0.08)
    ylabel("Frequency", fontsize="xx-large")
    xlabel("Average I₂→E Assembly Coupling Strength (pF)", fontsize="xx-large")
    nbins = 60
    alpha = 0.2 # opacity
    ln_width = 3.5
    flt_window = 2.5
    hcolors = matplotlib.cm.Reds(56:round(Int, 200/7):256)
    pcolors = copy(hcolors)
    hcolors[:, 4] .= alpha
    # ii-16.2
    n, bins, patches = hist(x16, nbins, density=true, fc=hcolors[1, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[1, :], label=string("I₁→I₂=16.2, Q=", round(q16, digits=2), ", (μ=", round(mean(x16), digits=1), ", σ=", round(std(x16), digits=1), ")"), linewidth=ln_width)
    # ii-20
    n, bins, patches = hist(x20, nbins, density=true, fc=hcolors[2, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[2, :], label=string("I₁→I₂=20, Q=", round(q20, digits=2), ", (μ=", round(mean(x20), digits=1), ", σ=", round(std(x20), digits=1), ")"), linewidth=ln_width)
    # ii-28
    n, bins, patches = hist(x28, nbins, density=true, fc=hcolors[3, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[3, :], label=string("I₁→I₂=28, Q=", round(q28, digits=2), ", (μ=", round(mean(x28), digits=1), ", σ=", round(std(x28), digits=1), ")"), linewidth=ln_width)
    # ii-30
    n, bins, patches = hist(xOrig, nbins, density=true, fc=hcolors[4, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[4, :], label=string("*I₁→I₂=30, Q=", round(qOrig, digits=2), ", (μ=", round(mean(xOrig), digits=1), ", σ=", round(std(xOrig), digits=1), ")"), linewidth=ln_width)
    # ii-32
    n, bins, patches = hist(x32, nbins, density=true, fc=hcolors[5, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[5, :], label=string("I₁→I₂=32, Q=", round(q32, digits=2), ", (μ=", round(mean(x32), digits=1), ", σ=", round(std(x32), digits=1), ")"), linewidth=ln_width)
    # τ200
    n, bins, patches = hist(x40, nbins, density=true, fc=hcolors[6, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[6, :], label=string("I₁→I₂=40, Q=", round(q40, digits=2), ", (μ=", round(mean(x40), digits=1), ", σ=", round(std(x40), digits=1), ")"), linewidth=ln_width)
    # τ500
    n, bins, patches = hist(x50, nbins, density=true, fc=hcolors[7, :])
    gaus_kernel = ImageFiltering.Kernel.gaussian((flt_window,))
    data = imfilter(n, gaus_kernel)
    plot(bins, [data; 0], color=pcolors[7, :], label=string("I₁→I₂=50, Q=", round(q50, digits=2), ", (μ=", round(mean(x50), digits=1), ", σ=", round(std(x50), digits=1), ")"), linewidth=ln_width)

    ax = gca()
    lgd = legend(loc="upper right", fontsize="x-large")
    ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
    savefig(string("weight_dist_ii.png"), bbox_inches="tight")
    PyPlot.clf()
end

function getSequenceActivations(assemAct, seq_length, Npop; time_limit=100, overlap_limit=20)
    # Check for successful activation of whole sequences
    sequenceAct = zeros(1, 3)
    for ipop = 1:seq_length:Npop
        active_times = assemAct[assemAct[:, 1] .== ipop, 2:3]
        for (ind, tt) in enumerate(active_times[:, 1])
            end_time = successSequence(assemAct, active_times[ind, 2], seq_length-1, time_limit, overlap_limit, ipop+1)
            if end_time > 0
                sequenceAct = [sequenceAct; [ipop tt end_time]]
            end
        end
    end
    return sequenceAct[2:end, :]
end

function successSequence(assemAct, endActTime, seq_length, time_timit, overlap_limit, ipop)
    if seq_length == 0
        return endActTime
    else
        active_times = assemAct[assemAct[:, 1] .== ipop, 2:3]
        if isempty(active_times)
            return 0
        else
            for (ind, tt) in enumerate(active_times[:, 1])
                if ((tt - endActTime) <= time_limit) && ((tt - endActTime) >= -overlap_limit)
                    return successSequence(assemAct, active_times[ind, 2], seq_length-1, time_limit, overlap_limit, ipop+1)
                end
            end
            return 0
        end
    end
end

function getAssemblyActivations(eAssemfRates, Npop; threshold=10, fRate_window=2)
    gaus_kernel = ImageFiltering.Kernel.gaussian((fRate_window,))
    sucAct = zeros(size(eAssemfRates))
    for ipop = 1:Npop
        sucAct[ipop, :] = (imfilter(eAssemfRates[ipop, :], gaus_kernel) .> fr_threshold)
    end
    success_flag = zeros(Npop)
    assemAct = zeros(1, 3)
    for ipop = 1:Npop
        for tt = 1:size(sucAct)[2]
            if sucAct[ipop, tt] == 1 && success_flag[ipop] == 0
                success_flag[ipop] = tt
            end
            if sucAct[ipop, tt] == 0 && success_flag[ipop] != 0
                assemAct = [assemAct; [ipop success_flag[ipop] tt]]
                success_flag[ipop] = 0
            end
        end
    end
    return assemAct[2:end, :]
end

function getAssemblySpikes(Nsteps, Npop, popmembers, ipopmembers, times, dt)
    # Returns the total spike count through time (timesteps) for each excitatory and inhibitory assembly
    eAssemSpikes = zeros(Npop, Nsteps)  # Contains total spikes for each e-assembly at each timestep
    iAssemSpikes = zeros(Npop, Nsteps)  # Contains total spikes for each i-assembly at each timestep
    for ipop = 1:Npop
        members = filter(i->(i>0), popmembers[ipop, :])
        for member in members
            for tt in times[member, :]
                (tt == 0) && continue
                t = round(Int, tt/dt)
                eAssemSpikes[ipop, t] += 1
            end
        end
        imembers = filter(i->(i>0), ipopmembers[ipop, :])
        for imember in imembers
            for tt in times[imember, :]
                (tt == 0) && continue
                t = round(Int, tt/dt)
                iAssemSpikes[ipop, t] += 1
            end
        end

    end
    return eAssemSpikes, iAssemSpikes
end
