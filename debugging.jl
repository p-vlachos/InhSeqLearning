sim_name = string("cross_corr_spontaneous.h5")
sim_savedpath = "./analysis_data/"

fid = h5open(joinpath(sim_savedpath, sim_name), "r")
crossEE = read(fid["data"]["crossEE"])
crossEI = read(fid["data"]["crossEI"])
crossIE = read(fid["data"]["crossIE"])
close(fid)


# Average over all simulations
avg_cross = mean(crossEE, dims=4)
avg_cross = crossEE[:, :, :, 7]
# Average over all 1st sequences
cross0 = mean(avg_cross[:, 1, 1, 1], dims=2)
cross0 += mean(avg_cross[:, 5, 5, 1], dims=2)
cross0 += mean(avg_cross[:, 9, 9, 1], dims=2)
cross0 += mean(avg_cross[:, 13, 13, 1], dims=2)
cross0 += mean(avg_cross[:, 17, 17, 1], dims=2)

cross1 = mean(avg_cross[:, 1, 2, 1], dims=2)
cross1 += mean(avg_cross[:, 5, 6, 1], dims=2)
cross1 += mean(avg_cross[:, 9, 10, 1], dims=2)
cross1 += mean(avg_cross[:, 13, 14, 1], dims=2)
cross1 += mean(avg_cross[:, 17, 18, 1], dims=2)

cross2 = mean(avg_cross[:, 1, 3, 1], dims=2)
cross2 += mean(avg_cross[:, 5, 7, 1], dims=2)
cross2 += mean(avg_cross[:, 9, 11, 1], dims=2)
cross2 += mean(avg_cross[:, 13, 15, 1], dims=2)
cross2 += mean(avg_cross[:, 17, 19, 1], dims=2)

cross3 = mean(avg_cross[:, 1, 4, 1], dims=2)
cross3 += mean(avg_cross[:, 5, 8, 1], dims=2)
cross3 += mean(avg_cross[:, 9, 12, 1], dims=2)
cross3 += mean(avg_cross[:, 13, 16, 1], dims=2)
cross3 += mean(avg_cross[:, 17, 20, 1], dims=2)



legendlabelsize = 28
labelsize = 36
ticklabelsize = 36

baseline = 17
fig = Figure(resolution = (1080, 720))
ax_single = Axis(fig[1, 1], 
            xlabel=L"\text{lag (ms)}", xlabelsize=labelsize,
            ylabel=L"\text{cross-correlation}", ylabelsize=labelsize,
            xticks=([-200., -100., 0., 100., 200.], [L"-200", L"-100", L"0", L"100", L"200"]), xticklabelsize=ticklabelsize,
            yticks=([0., .5, 1.], [L"0", L"0.5", L"1"]), yticklabelsize=ticklabelsize,
            limits=(-280, 280, -0.15, 1.1))
#1
for pop = 20:-1:5
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end
for pop = 1:4
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Blues[5+pop], linewidth=linewidth/2, label=L"\text{members}")
end

#5
for pop = 20:-1:9
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end
for pop = 5:9
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Blues[pop], linewidth=linewidth/2, label=L"\text{members}")
end
for pop = 4:-1:1
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end

#9
for pop = 20:-1:13
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end
for pop = 9:12
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Blues[pop-4], linewidth=linewidth/2, label=L"\text{members}")
end
for pop = 8:-1:1
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end

#13
for pop = 20:-1:17
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end
for pop = 13:16
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Blues[pop-8], linewidth=linewidth/2, label=L"\text{members}")
end
for pop = 12:-1:1
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end


#17
for pop = 17:20
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Blues[pop-12], linewidth=linewidth/2, label=L"\text{members}")
end
for pop = 16:-1:1
    lines!(ax_single, -400:400, avg_cross[:, baseline, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{other}")
end

# lines!(ax_single, -400:400, avg_cross[:, 1, 1], color=ColorSchemes.Blues[6], linewidth=linewidth, label=L"E1-E1")
# lines!(ax_single, -400:400, avg_cross[:, 1, 2], color=ColorSchemes.Blues[7], linewidth=linewidth, label=L"E1-E2")
# lines!(ax_single, -400:400, avg_cross[:, 1, 3], color=ColorSchemes.Blues[8], linewidth=linewidth, label=L"E1-E3")
# lines!(ax_single, -400:400, avg_cross[:, 1, 4], color=ColorSchemes.Blues[9], linewidth=linewidth, label=L"E1-E4")
axislegend(ax_single, position=(0.95, 0.95), labelsize=legendlabelsize, merge=true)
fig











using CairoMakie
using ColorSchemes


cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))

cl_inh = ColorScheme(range(colorant"gray80", colorant"firebrick", length=100))
cl_exc = ColorScheme(range(colorant"gray80", colorant"dodgerblue4", length=100))

inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))
excitation_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_exc, LinRange(0, 1, 100)))


sim_name = string("network_7_stimulation.h5")
sim_savedpath = "./networks_trained_stimulation/"
output_dir = "./output_analysis/"

fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
weightsEE = read(fid["data"]["weightsEE"])
weightsEI = read(fid["data"]["weightsEI"])
weightsIE = read(fid["data"]["weightsIE"])
times = read(fid["data"]["times"])
close(fid)



# plotNetworkActivity(times, popmembers, ipopmembers; interval=5_000:15_000, name="_testin")
plotWeightsEE(meanEE, name="_testinEE")
plotWeightsIE(mean(weightsIE[ipopmembers[:, :] .- 4000, :, 1000], dims=1)[1, :, :], name="_testinIE")
plotWeightsEI(mean(weightsEI[ipopmembers[:, :] .- 4750, :, 1000], dims=1)[1, :, :], name="_testinEI")

size(weightsEI[ipopmembers[:, :] .- 4750, :, 1000])

ipopmembers = findI2populations(weights, 20, popmembers, iipop_len=27)

meanEE = zeros(20, 20)
for ipop = 1:20
	for iipop = 1:20
		meanEE[ipop, iipop] = sum(weights[filter(i->i>0, popmembers[:, ipop]), filter(i->i>0, popmembers[:, iipop])]) / count(i->i>0, weights[filter(i->i>0, popmembers[:, ipop]), filter(i->i>0, popmembers[:, iipop])]) 
	end
end

meanEI = zeros(20, 20)
for ipop = 1:20
	for iipop = 1:20
		meanEI[ipop, iipop] = mean(weightsEI[ipopmembers[:, ipop] .- 4750, iipop, 1000])
	end
end

meanIE = zeros(20, 20)
for ipop = 1:20
	for iipop = 1:20
		meanIE[ipop, iipop] = mean(weightsIE[ipopmembers[:, ipop] .- 4000, iipop, 1000])
	end
end

heatmap(meanEI, colormap=excitation_cs)
heatmap(meanIE, colormap=inhibition_cs)


simulation = 2
sim_name = string("network", simulation, "_.h5")
sim_savedpath = "./networks_trained/"
sim_savedpath = "./networks_trained_spontaneous/"
sim_savedpath = "./networks_trained_stimulation/"
for sim = simulation:simulation

	if sim_savedpath == "./networks_trained/"
		sim_name = string("network_", sim,".h5")
	elseif sim_savedpath == "./networks_trained_spontaneous/"
		sim_name = string("network_", sim,"_spontaneous.h5")
	else
		sim_name = string("network_", sim,"_stimulation.h5")
	end

	fid = h5open(joinpath(sim_savedpath, sim_name), "r")
	popmembers = read(fid["data"]["popmembers"])
	weights = read(fid["data"]["weights"])
	weightsEE = read(fid["data"]["weightsEE"])
	weightsEI = read(fid["data"]["weightsEI"])
	weightsIE = read(fid["data"]["weightsIE"])
	times = read(fid["data"]["times"])
	close(fid)

	# ipopmembers = findI2populations(weights, 20, popmembers, iipop_len=27)
	ipopmembers = findI2populations(weights, 16, popmembers, iipop_len=27)
	output_dir = string("./output_analysis/simulation_", sim, "/") 
	# plotWeightsEE(meanEE, name="_testinEE", output_dir=output_dir)
	if sim_savedpath == "./networks_trained/"
		plotWeightsEE(weightsEE[:, :, 1000], name="_testinEE", output_dir=output_dir)
		plotWeightsIE(mean(weightsIE[ipopmembers[:, :] .- 4000, :, 1000], dims=1)[1, :, :], name="_testinIE", output_dir=output_dir)
		plotWeightsEI(mean(weightsEI[ipopmembers[:, :] .- 4750, :, 1000], dims=1)[1, :, :], name="_testinEI", output_dir=output_dir)
	elseif sim_savedpath == "./networks_trained_spontaneous/"
		plotNetworkActivity(times, popmembers, ipopmembers; interval=5_000:10_000, name=string("_testinActivitySpontaneous_", sim))
	else
		plotNetworkActivity(times, popmembers, ipopmembers; interval=5_000:10_000, name=string("_testinActivityStimulation_", sim))
	end
end






sim_name = string("network_1_spontaneous.h5")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
weightsEE = read(fid["data"]["weightsEE"])
weightsEI = read(fid["data"]["weightsEI"])
weightsIE = read(fid["data"]["weightsIE"])
times = read(fid["data"]["times"])
close(fid)





###########################################################################################
##########                              RATE MODEL                    			 ##########
###########################################################################################
using Distributions
using ProgressBars
using Random


T = 1000

w, x = rate_simnew(T)

fig = Figure(resolution=(720, 480))
ax = Axis(fig[1, 1], xlabel="time (a.u.)", ylabel="firing rate (a.u.)")

for ipop = 1:4
	lines!(ax, x[ipop], color="blue")
end
fig

for ipop = 6:9
	lines!(ax, x[ipop], color="red")
end

fig
lines!(ax, x[5], color="green")
fig
lines!(ax, x[10], color="black")
fig
# for ipop = 1:4
# 	lines!(ax, Φ.(x[ipop]))
# end

fig
fig














###########################################################################################
##########                          Cross-Correlation                 			 ##########
###########################################################################################
# using StatsBase
using DSP
using CairoMakie

function gen_gaussian(matrix, center, sigma)
	# Generate a 2D gaussian at "center" with specified "sigma"
	x_c, y_c = center
	for x in 1:size(matrix, 1)
		for y in 1:size(matrix, 2)
			matrix[x, y] += exp(-((x - x_c)^2 + (y - y_c)^2) / (2*sigma^2))
		end
	end
	return matrix
end


x = zeros(1024, 1024)
x = gen_gaussian(x, (300, 400), 80.)
x = gen_gaussian(x, (750, 800), 90.)
x = 1 .- x
heatmap(x, colormap="Greys")

y = zeros(1024, 1024)
y = gen_gaussian(y, (200, 200), 80.)
y = gen_gaussian(y, (800, 900), 90.)
y = 1 .- y
heatmap(y, colormap="Greys")


function comp_crosscor(x, y)
	ksize = 32
	overlap = 8
	# z = zeros()
	for i = 1:(ksize-overlap):(1024-ksize)
		for j = 1:(ksize-overlap):(1024-ksize)
			crosscor(x[i:i+ksize, j:j+ksize], y[i:i+ksize, j:j+ksize])
			# sum([x[i:i+ksize, j:j+ksize], y[i:i+ksize, j:j+ksize]])
		end
	end
end

@time comp_crosscor(x, y)


[vec(x[i:i+ksize, j:j+ksize]); vec(y[i:i+ksize, j:j+ksize])]

z = DSP.xcorr(x[i:i+ksize, j:j+ksize], y[i:i+ksize, j:j+ksize])


x[i:i+ksize, j:j+ksize]

z = crosscor(vec(x[i:i+ksize, j:j+ksize]), vec(y[i:i+ksize, j:j+ksize]))

heatmap(z[3, :, :])


a = rand(10)
b = rand(10)

c = crosscor(a, b)

lines(z)


i = j = 1
ksize=32
overlap=8
size(1:(ksize-overlap):(1024-ksize))


heatmap(x, colormap="Greys")
heatmap(y, colormap="Greys")



































# ________________________________________________________________________
T = 8_000
stim = [1. 1001 1030 8; 1 1201 1230 8; 1 1401 1430 8; 1 1601 1630 8; 1 1801 1830 8;
		5 2001 2030 8; 5 2201 2230 8; 5 2401 2430 8; 5 2601 2630 8; 5 2801 2830 8;
		9 3001 3030 8; 9 3201 3230 8; 9 3401 3430 8; 9 3601 3630 8; 9 3801 3830 8;
		13 4001 4030 8; 13 4201 4230 8; 13 4401 4430 8; 13 4601 4630 8; 13 4801 4830 8;
		17 5001 5030 8; 17 5201 5230 8; 17 5401 5430 8; 17 5601 5630 8; 17 5801 5830 8]

sim_name = string("network_1.h5")
sim_savepath = "./networks/"
output_dir = "./output/"

fid = h5open(joinpath(sim_savepath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights_old = read(fid["data"]["weights"])
close(fid)
times, ns, Ne, Ncells, T, new_weights = sim(stim, weights_old, popmembers, T)


ipopmembers = findI2populations(weights_old, 20, popmembers, iipop_len=25)

i_structs = zeros(20, 500)

for i = 1:20
	members = filter(i->i>0, popmembers[i, :])
	for j = 4501:5000
		# i_structs[i, 1:length(weights_old[j, members])] .= weights_old[j, members]
		# i_structs[i, 1:length(members)] .= weights_old[j, members]
		for mem in members
			i_structs[i, j-4500] += weights_old[j, mem]
		end
	end
end

i_count = zeros(20, 500)

for i = 1:20
	members = filter(i->i>0, popmembers[i, :])
	for j = 4501:5000
		for mem in members
			if weights_old[j, mem] > 0
				i_count[i, j-4500] += 1
			end
		end
	end
end


for i = 1:20
	perm = sortperm(i_structs[2, :], rev=true)
	if i==1
		Plots.plot(i_structs[i, perm])
		# Plots.plot(sort(i_structs[i, :], rev=true))
		# Plots.plot(i_structs[i, :])
	else
		Plots.plot!(i_structs[i, perm])
		# Plots.plot!(sort(i_structs[i, :], rev=true))
		# Plots.plot!(i_structs[i, :])
	end
end

Plots.plot!(legend=false)








###############################################################################
ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=25)

totalI = zeros(20)
for ipop = 1:20
	for (imem, mem) in enumerate(popmembers[ipop, :])
		if mem == 0
			continue
		end
		for i = 4500:5000
			totalI[ipop] += w8s[mem, i]
		end
	end
end

###############
using Clustering

# probably I have to flip
i_rates = convolveSpikes(times[(Ncells-Ni2+1):Ncells, :], interval=1:T)
# i_rates = convolveSpikes(x, interval=1:T)

R = kmeans(i_rates', 21)
R = hclust(i_rates')

@assert nclusters(R) == 20 # verify the number of clusters

a = assignments(R) # get the assignments of points to clusters
c = counts(R) # get the cluster sizes
M = R.centers # get the cluster centers

# plot with the point color mapped to the assigned cluster index
scatter(1:500, a)



ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=25)
imshow(ipopmembers);colorbar()
overlap = zeros(20, 20)
for i = 1:20
	for j = 1:20
		if j == i
			continue
		end
		for zi = 1:25
			for zj = 1:25
				if ipopmembers[i, zi] == ipopmembers[j, zj]
					overlap[i, j] += 1
				end
			end
		end
	end
end

misin = 0
for i in eachindex(ipopmembers)
	if ipopmembers[i] < 4501 || ipopmembers[i] > 5000
		misin += 1
	end
end

# --- Old Stuff --- 

similar = zeros(20, 20)

for pop = 1:20
	for imember in ipops[pop, :]
		for pop2 = 1:20
			(pop2 == pop) && (continue)
			for imember2 in ipops[pop2, :]
				if imember == imember2
					similar[pop, pop2] += 1
				end
			end
		end
	end
end


similar2 = zeros(20, 20)

for pop = 1:20
	for imember in ipops2[pop, :]
		for pop2 = 1:20
			(pop2 == pop) && (continue)
			for imember2 in ipops2[pop2, :]
				if imember == imember2
					similar2[pop, pop2] += 1
				end
			end
		end
	end
end


similar3 = zeros(20, 20)

for pop = 1:20
	for imember in ipops3[pop, :]
		for pop2 = 1:20
			(pop2 == pop) && (continue)
			for imember2 in ipops3[pop2, :]
				if imember == imember2
					similar3[pop, pop2] += 1
				end
			end
		end
	end
end


similare1 = zeros(20, 20)

for pop = 1:20
	for imember in pops[pop, :]
		(imember == 0) && (continue)
		for pop2 = 1:20
			(pop2 == pop) && (continue)
			for imember2 in pops[pop2, :]
				(imember == 0) && (continue)
				if imember == imember2
					similare1[pop, pop2] += 1
				end
			end
		end
	end
end


similare2 = zeros(20, 20)

for pop = 1:20
	for imember in pops2[pop, :]
		(imember == 0) && (continue)
		for pop2 = 1:20
			(pop2 == pop) && (continue)
			for imember2 in pops2[pop2, :]
				(imember == 0) && (continue)
				if imember == imember2
					similare2[pop, pop2] += 1
				end
			end
		end
	end
end

similare3 = zeros(20, 20)

for pop = 1:20
	for imember in pops3[pop, :]
		(imember == 0) && (continue)
		for pop2 = 1:20
			(pop2 == pop) && (continue)
			for imember2 in pops3[pop2, :]
				(imember == 0) && (continue)
				if imember == imember2
					similare3[pop, pop2] += 1
				end
			end
		end
	end
end






similari = zeros(20, 1000)

for pop = 1:20
	for imember in pops2[pop, :]
		(imember == 0) && (continue)
		for imember2 = 4001:5000
			(w8s2[imember, imember2] == 0) && (continue)
			similari[pop, imember2-4000] += 1
		end
	end
end








epopsums2 = zeros(20)
for ipop = 1:20
	emems = filter(i->(i>0), pops2[ipop, :])
	epopsums2[ipop] = sum(w8s2[4001:4500, emems])
	# epopsums2[ipop] = sum(w8s2[4501:5000, emems])
end


avg_E3 = zeros(20)
for ipop = 1:20
	emems = filter(i->(i>0), pops3[ipop, :])
	avg_E3[ipop] = sum(w8s3[emems, emems])
end


output_dir = "./output/testing/"
sim_outDir = string(output_dir, "/simtest/")
plot_eeWeights(w8s3, pops3, 20, sim_outDir)
plot_eiWeights(w8s, pops, ipops, 20, sim_outDir)

##### =========  OLD STUFF  ===========
# overlap = calc_overlap(popmembers, Ne, Npop)
ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=25)
eAssemSpikes, iAssemSpikes = getAssemblySpikes(Nsteps, Npop, popmembers, ipopmembers, times, dt)
eAssemfRates, iAssemfRates = getAssemblyfRates(Npop, T, eAssemSpikes, iAssemSpikes)
assemAct = getAssemblyActivations(eAssemfRates, Npop, threshold=fr_threshold, fRate_window=fRate_window)
sequenceAct = getSequenceActivations(assemAct, seq_length, Npop, time_limit=time_limit, overlap_limit=overlap_limit)

plot_eeWeights(popmembers, Npop)
plot_eiWeights(popmembers, ipopmembers, Npop)
plot_ieWeights(popmembers, ipopmembers, Npop)
plot_eiWeights_zoom(popmembers, ipopmembers, Npop, seq_num=1, seq_len=seq_length)
plot_eiWeights_zoom2(popmembers, ipopmembers, Npop, seq_num=3, seq_len=seq_length)


# barplot_eiRecWeights(popmembers, ipopmembers, Npop)
# barplot_ieRecWeights(popmembers, ipopmembers, Npop)
# barplot_eiRecWeights(popmembers2, ipopmembers, Npop)
# compare_eRecWeights(popmembers, popmembers2, Npop)
# compare_ieRecWeights(popmembers, ipopmembers, Npop, popmembers2, ipopmembers2)
# compare_ieRecWeights2(popmembers, ipopmembers, Npop, popmembers2, ipopmembers2)
# compare_eiRecWeights(popmembers, ipopmembers, Npop, popmembers2, ipopmembers2)
# plot_AssemfRate(Npop, T, eAssemfRates, iAssemfRates, seq_len=seq_length, fRate_window=fRate_window)

# plotSingleTotalInput(Ncells, times, T)

ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=25)
eRates = convolveSpikes(times, interval=1:50000, tau=10)
eAssemfRates, iAssemfRates, restfRates, inh1fRates = getAssemblySpikes1(eRates, popmembers, ipopmembers)

excfRates = sum(eAssemfRates, dims=1)' ./ 20
inh2fRates = sum(iAssemfRates, dims=1)' ./ 20
allifRates = sum(eRates[4000:5000, :], dims=1)' ./ 1000


plotGroupfRates2(vec(excfRates), vec(inh1fRates), vec(inh2fRates), vec(restfRates))
print()


plotGroupfRates(excfRates, inh1fRates, inh2fRates, mint=1000, maxt6=6000)
plotGroupfRates(excfRates, allifRates, mint=1000, maxt=6000)
plot_AssemfRate(Npop, T, eAssemfRates, iAssemfRates)

eAssemAct = findAssemblyActivations(eAssemfRates)
iAssemAct = findAssemblyActivations(iAssemfRates)

succeSeqAct = findSeqActivations(eAssemAct)
succiSeqAct = findSeqActivations(iAssemAct)

succeSeqCount = numSeqActivations(succeSeqAct)
succiSeqCount = numSeqActivations(succiSeqAct)

asseAct, correCoef = seqCorrelationCoef(eAssemAct)
assiAct, corriCoef = seqCorrelationCoef(iAssemAct)


assieAct, corrIECoef = seqEICorrelationCoef(eAssemAct, iAssemAct)

# eRates500 = copy(eRates)
# w8s500 = copy(new_weights)
# popmembers500 = copy(popmembers)
# ipopmembers500 = copy(ipopmembers)

dataEI, dataIE, summed, NsuccAct, NfailAct, NassemAct = plot_seqCorrelationCoef(eAssemAct, popmembers, ipopmembers)

corrIECoef, corrEICoef = plot_seqCorrelationCoef(eAssemAct, popmembers, ipopmembers)
# NOTE: HERE a value of 1 means that in 100% of cases were sequential dynamics were achieved, inhibition was responsible for it

qValue, variance = calculateQ(eAssemAct)



eAssemfRates, iAssemfRates, restfRates, inh1fRates = getAssemblySpikes1(eRatesHigh, popmembersHigh, ipopmembersHigh)
eAssemfRates, iAssemfRates, restfRates, inh1fRates = getAssemblySpikes1(eRatesLow, popmembersLow, ipopmembersLow)
# y50 = [mean(eAssemfRates) mean(iAssemfRates) mean(restfRates) mean(inh1fRates)]

#############################################################################################

if true	# TRAINED
	pygui(false)
	Npop = size(popmembers, 1)
	Nmaxmembers = size(popmembers, 2)
	filter_min = 16000
	filter_max = 21000
	ipopsize = 25
	labelsize = 15.5
	axessize = "x-large"
	smallaxessize = "large"
	ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=ipopsize)
	restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
	println("Creating plot...")
	fig, axs = subplots(2, 1, figsize=(7, 8), gridspec_kw=Dict("height_ratios"=>[1, 4]), sharex=true)
	axs[2].set_xlim(filter_min, filter_max)
	axs[2].set_ylim(0, sum(popmembers .> 0)+500+sum(ipopmembers .> 0)+length(restcells))
	lineWidth = .3
	ylabel("Sequences (Neurons)", fontsize=labelsize, labelpad=6)
	xlabel("Simulation Time (ms)", fontsize=labelsize)
	# Plot raster with the order of rows determined by population membership
	rowcount = 0
	ytickLoc = zeros(22)
	for pp = 1:Npop
		print("\rpopulation ", pp)
		for cc = 1:Nmaxmembers
			if popmembers[pp, cc] < 1
				break
			end
			global rowcount += 1
			ind = popmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			axs[2].scatter(vals, y, s = .1, c="tab:blue", marker="o", linewidths=lineWidth, label="Eₘ")
		end
		ytickLoc[pp] = rowcount
	end
	locEm = rowcount / 2
	# Find and plot raster of cells not belonging to any population
	for cc in restcells
		global rowcount += 1
		global vals = times[cc, 1:ns[cc]]
		global vals = filter(i->filter_max>i>filter_min, vals)
		global y = rowcount * ones(length(vals))
		axs[2].scatter(vals, y, s = .1, c="tab:gray", marker="o", linewidths=0.3, label="Eᵣ")
	end
	ytickLoc[Npop+1] = rowcount
	locEn = (2 * locEm) + ((rowcount - 2 * locEm) / 2)
	for cc = 4001:4500
	# for cc = 4001:5000
		global rowcount += 1
		global vals = times[cc, 1:ns[cc]]
		global vals = filter(i->filter_max>i>filter_min, vals)
		global y = rowcount * ones(length(vals))
		axs[2].scatter(vals, y, s = .1, c="brown", marker="o", linewidths=lineWidth, label="I₁")
	end
	ytickLoc[Npop+2] = rowcount
	locI1 = (2 * (locEn-locEm)) + ((rowcount - (2 * (locEn-locEm))) / 2)
	# Plot raster for interneurons (I2)
	for pp = 1:Npop
		print("\ri-population ", pp)
		for cc = 1:ipopsize
			global rowcount += 1
			ind = ipopmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			axs[2].scatter(vals, y, s = .1, c="red", marker="o", linewidths=lineWidth, label="I₂")
		end
	end
	locI2 = (2 * (locI1 - locEn + locEm)) + ((rowcount - (2 * (locI1 - locEn + locEm))) / 2)

	minor_ticks = matplotlib.ticker.FixedLocator(ytickLoc)
	assem_ticks = zeros(round(Int, Npop/seq_length))
	iassem_ticks = zeros(round(Int, Npop/seq_length))
	cc = 0
	for tick = seq_length:seq_length:Npop
		cc += 1
		(cc == 1) ? assem_ticks[cc] = ytickLoc[tick] / 2 : assem_ticks[cc] = (ytickLoc[tick-seq_length] + (ytickLoc[tick] - ytickLoc[tick-seq_length])/ 2)
	end

	major_ticks = matplotlib.ticker.FixedLocator(assem_ticks)
	axs[2].yaxis.set_major_locator(major_ticks)
	axs[2].set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize=axessize, ha="center", fontfamily="monospace")
	axs[2].yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false)
	axs[2].yaxis.set_minor_locator(minor_ticks)
	axs[2].yaxis.set_tick_params(which="minor", labelright=false, labelleft=false, right=false, left=true)
	for pp = 1:Npop
		if mod(pp,seq_length) == 0
			axs[2].axhline([ytickLoc[pp]], linestyle="-", linewidth=.5, color="black")
			axs[2].axhline([ytickLoc[end]+(ipopsize * pp)], linestyle="-", linewidth=.3, color="black")
		else
			axs[2].axhline([ytickLoc[pp]], linestyle="-", linewidth=.15, color="black")
		end
	end
	axs[2].axhline([ytickLoc[Npop+1]], linestyle="-", linewidth=.5, color="black");axs[2].axhline([ytickLoc[Npop+2]], linestyle="-", linewidth=.5, color="black")
	
	ax2 = axs[2].secondary_yaxis("left")
	major_ticks = matplotlib.ticker.FixedLocator([(locI2-2*(ipopsize*seq_length)), (locI2-(ipopsize*seq_length)), locI2, (locI2+(ipopsize*seq_length)), (locI2+(2*ipopsize*seq_length))])
	ax2.yaxis.set_major_locator(major_ticks)
	ax2.set_yticklabels(["A", "B  ", "C", "D  ", "E"], minor=false, fontsize=smallaxessize, ha="center", fontfamily="monospace")
	ax2.yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false)

	xticks = matplotlib.ticker.FixedLocator(filter_min:1000:filter_max)
	axs[2].xaxis.set_major_locator(xticks)
	axs[2].set_xticklabels(0:1000:(filter_max-filter_min), minor=false, fontsize=axessize)
	
	# ___ --- PLOT THE FIRING RATES --- ___
	print("\rPlotting firing rates")

	eRates = convolveSpikes(times, interval=filter_min:filter_max, tau=10)
	eAssemfRates, iAssemfRates, restfRates, inh1fRates = getAssemblySpikes1(eRates, popmembers, ipopmembers)
	excfRates = sum(eAssemfRates, dims=1)' ./ 20
	inh2fRates = sum(iAssemfRates, dims=1)' ./ 20

    lineWidth = 1
    axs[1].set_ylabel("Firing Rate (Hz)", fontsize=labelsize)

    axs[1].plot([zeros(filter_min); vec(excfRates)], color="tab:blue", label="Excitatory (members)", lw=lineWidth)
    axs[1].plot([zeros(filter_min); vec(restfRates)], color="tab:gray", label="Excitatory (non-members)", lw=lineWidth)
    axs[1].plot([zeros(filter_min); vec(inh1fRates)], color="brown", label="Inhibitory₁", lw=lineWidth)
    axs[1].plot([zeros(filter_min); vec(inh2fRates)], color="red", label="Inhibitory₂", lw=lineWidth)

    axs[1].xaxis.set_tick_params(labeltop=false, labelbottom=false, top=false, bottom=true, labelsize=smallaxessize)
    axs[1].yaxis.set_tick_params(labelright=false, labelleft=true, right=false, left=true, labelsize=smallaxessize)
	yticks = matplotlib.ticker.FixedLocator(0:4)
	axs[1].yaxis.set_major_locator(yticks)

	tight_layout()
	savefig(string("output_", (filter_max-filter_min), ".png"), dpi=150, bbox_inches="tight")
	print("\rDone creating plot\n")
	PyPlot.clf()
end

#############################################################################################

if true	# UNTRAINED
	pygui(false)
	Npop = size(popmembers, 1)
	Nmaxmembers = size(popmembers, 2)
	filter_min = 1000
	filter_max = 6000
	ipopsize = 25
	labelsize = 15.5
	axessize = "x-large"
	smallaxessize = "large"
	ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=ipopsize)
	restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
	println("Creating plot...")
	fig, axs = subplots(2, 1, figsize=(7, 8), gridspec_kw=Dict("height_ratios"=>[1, 4]), sharex=true)
	axs[2].set_xlim(filter_min, filter_max)
	axs[2].set_ylim(0, sum(popmembers .> 0)+500+sum(ipopmembers .> 0)+length(restcells))
	lineWidth = .3
	ylabel("Neurons", fontsize=labelsize, labelpad=6)
	xlabel("Simulation Time (ms)", fontsize=labelsize)
	# Plot raster with the order of rows determined by population membership
	rowcount = 0
	ytickLoc = zeros(22)
	for pp = 1:Npop
		print("\rpopulation ", pp)
		for cc = 1:Nmaxmembers
			if popmembers[pp, cc] < 1
				break
			end
			global rowcount += 1
			ind = popmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			axs[2].scatter(vals, y, s = .1, c="tab:blue", marker="o", linewidths=lineWidth, label="Eₘ")
		end
		ytickLoc[pp] = rowcount
	end
	locEm = rowcount / 2
	# Find and plot raster of cells not belonging to any population
	for cc in restcells
		global rowcount += 1
		global vals = times[cc, 1:ns[cc]]
		global vals = filter(i->filter_max>i>filter_min, vals)
		global y = rowcount * ones(length(vals))
		axs[2].scatter(vals, y, s = .1, c="tab:gray", marker="o", linewidths=0.3, label="Eᵣ")
	end
	ytickLoc[Npop+1] = rowcount
	locEn = (2 * locEm) + ((rowcount - 2 * locEm) / 2)
	for cc = 4001:4500
	# for cc = 4001:5000
		global rowcount += 1
		global vals = times[cc, 1:ns[cc]]
		global vals = filter(i->filter_max>i>filter_min, vals)
		global y = rowcount * ones(length(vals))
		axs[2].scatter(vals, y, s = .1, c="brown", marker="o", linewidths=lineWidth, label="I₁")
	end
	ytickLoc[Npop+2] = rowcount
	locI1 = (2 * (locEn-locEm)) + ((rowcount - (2 * (locEn-locEm))) / 2)
	# Plot raster for interneurons (I2)
	for pp = 1:Npop
		print("\ri-population ", pp)
		for cc = 1:ipopsize
			global rowcount += 1
			ind = ipopmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			axs[2].scatter(vals, y, s = .1, c="red", marker="o", linewidths=lineWidth, label="I₂")
		end
	end
	locI2 = (2 * (locI1 - locEn + locEm)) + ((rowcount - (2 * (locI1 - locEn + locEm))) / 2)

	minor_ticks = matplotlib.ticker.FixedLocator(ytickLoc)
	assem_ticks = zeros(round(Int, Npop/seq_length))
	iassem_ticks = zeros(round(Int, Npop/seq_length))
	cc = 0
	for tick = seq_length:seq_length:Npop
		cc += 1
		(cc == 1) ? assem_ticks[cc] = ytickLoc[tick] / 2 : assem_ticks[cc] = (ytickLoc[tick-seq_length] + (ytickLoc[tick] - ytickLoc[tick-seq_length])/ 2)
	end

	major_ticks = matplotlib.ticker.FixedLocator([locEm+550, locEn+250, locI1+150, locI2+100])
	axs[2].yaxis.set_major_locator(major_ticks)
	axs[2].set_yticklabels(["E-Members", "E-Rest", "Inh₁", "  Inh₂"], minor=false, fontsize=axessize, rotation=90, ha="right")
	axs[2].yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false, direction="out", pad=0.5)
	axs[2].yaxis.set_minor_locator(minor_ticks)
	axs[2].yaxis.set_tick_params(which="minor", labelright=false, labelleft=false, right=false, left=true)
	axs[2].axhline([ytickLoc[Npop]], linestyle="-", linewidth=.5, color="black")
	axs[2].axhline([ytickLoc[Npop+1]], linestyle="-", linewidth=.5, color="black")
	axs[2].axhline([ytickLoc[Npop+2]], linestyle="-", linewidth=.5, color="black")

	xticks = matplotlib.ticker.FixedLocator(filter_min:1000:filter_max)
	axs[2].xaxis.set_major_locator(xticks)
	axs[2].set_xticklabels(0:1000:(filter_max-filter_min), minor=false, fontsize=axessize)
	
	# ___ --- PLOT FIRING RATES --- ___
	print("\rPlotting firing rates")

	lineWidth = 1
	eRates = convolveSpikes(times, interval=filter_min:filter_max, tau=10)
	eAssemfRates, iAssemfRates, restfRates, inh1fRates = getAssemblySpikes1(eRates, popmembers, ipopmembers)
	excfRates = sum(eAssemfRates, dims=1)' ./ 20
	inh2fRates = sum(iAssemfRates, dims=1)' ./ 20

    axs[1].plot([zeros(filter_min); vec(excfRates)], color="tab:blue", label="Excitatory (members)", lw=lineWidth)
    axs[1].plot([zeros(filter_min); vec(restfRates)], color="tab:gray", label="Excitatory (non-members)", lw=lineWidth)
    axs[1].plot([zeros(filter_min); vec(inh1fRates)], color="brown", label="Inhibitory₁", lw=lineWidth)
    axs[1].plot([zeros(filter_min); vec(inh2fRates)], color="red", label="Inhibitory₂", lw=lineWidth)

    axs[1].set_ylabel("Firing Rate (Hz)", fontsize=labelsize)
    axs[1].xaxis.set_tick_params(labeltop=false, labelbottom=false, top=false, bottom=true, labelsize=smallaxessize)
    axs[1].yaxis.set_tick_params(labelright=false, labelleft=true, right=false, left=true, labelsize=smallaxessize)

	tight_layout()
	savefig(string("output_", (filter_max-filter_min), ".png"), dpi=150, bbox_inches="tight")
	print("\rDone creating plot\n")
	PyPlot.clf()
end

#########################################

if true
	pygui(false)
	Npop = size(popmembers, 1)
	Nmaxmembers = size(popmembers, 2)
	filter_min = 1000#000#1000#5000 #10300 #2200
	filter_max = 6000#6000 #10600 #2450
	ipopsize = 25
	ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=ipopsize)
	restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
	println("Creating plot...")
	figure(figsize=(7, 4))
	xlim(filter_min, filter_max)
	ylim(0, sum(popmembers .> 0)+500+sum(ipopmembers .> 0))#+length(restcells)
	lineWidth = .3
	ylabel("Sequences (Neurons)", fontsize=21, labelpad=6)
	xlabel("Simulation Time (ms)", fontsize=21)
	tight_layout()
	# Plot raster with the order of rows determined by population membership
	rowcount = 0
	ytickLoc = zeros(22)
	for pp = 1:Npop
		print("\rpopulation ", pp)
		for cc = 1:Nmaxmembers
			if popmembers[pp, cc] < 1
				break
			end
			global rowcount += 1
			ind = popmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			scatter(vals, y, s = .1, c="tab:blue", marker="o", linewidths=lineWidth, label="Eₘ")
		end
		ytickLoc[pp] = rowcount
	end
	locEm = rowcount / 2
	# Find and plot raster of cells not belonging to any population
	for cc in restcells
		global rowcount += 1
		global vals = times[cc, 1:ns[cc]]
		global vals = filter(i->filter_max>i>filter_min, vals)
		global y = rowcount * ones(length(vals))
		scatter(vals, y, s = .1, c="tab:gray", marker="o", linewidths=0.3, label="Eᵣ")
	end
	ytickLoc[Npop+1] = rowcount
	locEn = (2 * locEm) + ((rowcount - 2 * locEm) / 2)
	for cc = 4001:4500
	# for cc = 4001:5000
		global rowcount += 1
		global vals = times[cc, 1:ns[cc]]
		global vals = filter(i->filter_max>i>filter_min, vals)
		global y = rowcount * ones(length(vals))
		scatter(vals, y, s = .1, c="brown", marker="o", linewidths=lineWidth, label="I₁")
	end
	ytickLoc[Npop+2] = rowcount
	locI1 = (2 * (locEn-locEm)) + ((rowcount - (2 * (locEn-locEm))) / 2)
	# Plot raster for interneurons (I2)
	for pp = 1:Npop
		for cc = 1:ipopsize
			global rowcount += 1
			ind = ipopmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			scatter(vals, y, s = .1, c="red", marker="o", linewidths=lineWidth, label="I₂")
		end
	end
	locI2 = (2 * (locI1 - locEn + locEm)) + ((rowcount - (2 * (locI1 - locEn + locEm))) / 2)
	ax = gca()
	minor_ticks = matplotlib.ticker.FixedLocator(ytickLoc)
	assem_ticks = zeros(round(Int, Npop/seq_length))
	iassem_ticks = zeros(round(Int, Npop/seq_length))
	cc = 0
	for tick = seq_length:seq_length:Npop
		cc += 1
		(cc == 1) ? assem_ticks[cc] = ytickLoc[tick] / 2 : assem_ticks[cc] = (ytickLoc[tick-seq_length] + (ytickLoc[tick] - ytickLoc[tick-seq_length])/ 2)
	end

	# NOTE: USE THIS CODE FOR UNTRAINED NETWORKS
	# major_ticks = matplotlib.ticker.FixedLocator([locEm+550, locEn+250, locI1+150, locI2+100])
	# ax.yaxis.set_major_locator(major_ticks)
	# # ax.set_yticklabels(["Excitatory\n(members)", "Excitatory\n(non-members)", "Inhibitory₁", "Inhbitory₂"], minor=false, fontsize="x-small", rotation=45, ha="right")
	# # ax.set_yticklabels(["Eₘ", "Eᵣ", "I₁", "I₂"], minor=false, fontsize="x-small", rotation=45, ha="right")
	# ax.set_yticklabels(["E-Members", "E-Rest", "Inh₁", "Inh₂"], minor=false, fontsize="x-small", rotation=90, ha="right")
	# ax.yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false, direction="out", pad=0.5)
	# ax.yaxis.set_minor_locator(minor_ticks)
	# ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=false, right=false, left=true)
	# ax.axhline([ytickLoc[Npop]], linestyle="-", linewidth=.5, color="black")
	# ax.axhline([ytickLoc[Npop+1]], linestyle="-", linewidth=.5, color="black")
	# ax.axhline([ytickLoc[Npop+2]], linestyle="-", linewidth=.5, color="black")
	# #___________

	# NOTE: USE THIS CODE FOR TRAINED NETWORK
	major_ticks = matplotlib.ticker.FixedLocator(assem_ticks)
	ax.yaxis.set_major_locator(major_ticks)
	ax.set_yticklabels(["A", "B", "C", "D", "E"], minor=false, fontsize="xx-large", ha="center", fontfamily="monospace")
	ax.yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false)
	ax.yaxis.set_minor_locator(minor_ticks)
	ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=false, right=false, left=true)
	for pp = 1:Npop
		if mod(pp,seq_length) == 0
			ax.axhline([ytickLoc[pp]], linestyle="-", linewidth=.5, color="black")
			ax.axhline([ytickLoc[end]+(ipopsize * pp)], linestyle="-", linewidth=.3, color="black")
		else
			ax.axhline([ytickLoc[pp]], linestyle="-", linewidth=.15, color="black")
		end
	end
	ax.axhline([ytickLoc[Npop+1]], linestyle="-", linewidth=.5, color="black");ax.axhline([ytickLoc[Npop+2]], linestyle="-", linewidth=.5, color="black")
	ax2 = ax.secondary_yaxis("left")
	major_ticks = matplotlib.ticker.FixedLocator([(locI2-2*(ipopsize*seq_length)), (locI2-(ipopsize*seq_length)), locI2, (locI2+(ipopsize*seq_length)), (locI2+(2*ipopsize*seq_length))])
	ax2.yaxis.set_major_locator(major_ticks)
	ax2.set_yticklabels(["A", "B  ", "C", "D  ", "E"], minor=false, fontsize="large", ha="center", fontfamily="monospace")
	ax2.yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false)

	xticks = matplotlib.ticker.FixedLocator(filter_min:1000:filter_max)
	ax.xaxis.set_major_locator(xticks)
	ax.set_xticklabels(0:1000:(filter_max-filter_min), minor=false, fontsize="xx-large")

	savefig(string("output_", (filter_max-filter_min), ".png"), dpi=150)
	print("\rDone creating plot\n")
	PyPlot.clf()
end


plotGroupfRates(excfRates, inh1fRates, inh2fRates, mint=1000, maxt=6000)


# PLOT VOLTAGE THROUGH TIME
restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
emVolt = zeros(Nsteps); emInd = rand(filter(i->i>0, popmembers))
enVolt = zeros(Nsteps); enInd = rand(restcells)
i1Volt = zeros(Nsteps); i1Ind = rand(4001:4500)
i2Volt = zeros(Nsteps); i2Ind = rand(4501:5000)
for tt = 1:Nsteps
	emVolt[tt] = voltage_tracker[tt][emInd]
	enVolt[tt] = voltage_tracker[tt][enInd]
	i1Volt[tt] = voltage_tracker[tt][i1Ind]
	i2Volt[tt] = voltage_tracker[tt][i2Ind]
end

pygui(true)

for iseq = 1:Nseq
	plot(succeSeqAct[iseq,:])
end
plot(succeSeqAct[1,2000:4000])
plot(succeSeqAct[2,2000:4000])
plot(succeSeqAct[3,2000:4000])
plot(succeSeqAct[4,2000:4000])
plot(succeSeqAct[5,2000:4000])


for ipop = indexes
	plot(eAssemfRates[ipop,18600:18900])
end

for ipop = 17:17
	plot(iAssemfRates[ipop,18600:18900])
end

# for seq=1:5
if true
	pygui(false)
	filter_min = 8050
	filter_max = 8300
	seq = 5
	figure(figsize=(5, 5))
	xlim(filter_min, filter_max)
	indexes = 1:seq_length
	indexes = indexes .+ ((seq-1)*seq_length)
	ylim(0, 22)
	ylabel("E-Assembly\nFiring rate (Hz)", fontsize=30)
	# xlabel("Simulation Time (ms)", fontsize=30)
	ecolors = matplotlib.cm.Blues_r(1:round(Int,200/seq_length):200)
	tight_layout()
	for (ind, pp) in enumerate(indexes)
		plot(eAssemfRates[pp, :], color=ecolors[ind,:])
	end
	grid(b=true, which="major", axis="x", ls=":", c="silver")
	ax=gca()
	xticks = matplotlib.ticker.FixedLocator(filter_min:50:filter_max)
	ax.xaxis.set_major_locator(xticks)
	ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=false, top=false, bottom=true)
	ax.set_xticklabels(0:50:(filter_max-filter_min), minor=false, fontsize=23)
	ax.legend(["E-Assembly 17", "E-Assembly 18", "E-Assembly 19", "E-Assembly 20"], fontsize="xx-large", loc="upper right")
	y_ticks = matplotlib.ticker.FixedLocator(0:5:20)
	ax.yaxis.set_major_locator(y_ticks)
	yticks(fontsize=23)
	savefig(string("e_fr_", seq, ".png"),bbox_inches="tight", dpi=150)
	PyPlot.clf()

	xlim(filter_min, filter_max)
	ylim(0, 22)
	ylabel("I₂-Assembly\nFiring rate (Hz)", fontsize=30)
	xlabel("Simulation Time (ms)", fontsize=30)
	icolors = matplotlib.cm.Reds_r(1:round(Int,200/seq_length):200)
	tight_layout()
	for (ind, pp) in enumerate(indexes)
		plot(iAssemfRates[pp, :], color=icolors[ind,:])
	end
	grid(b=true, which="major", axis="x", ls=":", c="silver")
	ax=gca()
	xticks = matplotlib.ticker.FixedLocator(filter_min:50:filter_max)
	ax.xaxis.set_major_locator(xticks)
	# ax.set_xticklabels(0:50:(filter_max-filter_min), minor=false, fontsize=23)
	ax.set_xticklabels(["0", "", "100", "", "200", ""], minor=false, fontsize=23)
	ax.legend(["I₂-Assembly 17", "I₂-Assembly 18", "I₂-Assembly 19", "I₂-Assembly 20"], fontsize="xx-large", loc="upper right")
	y_ticks = matplotlib.ticker.FixedLocator(0:5:20)
	ax.yaxis.set_major_locator(y_ticks)
	yticks(fontsize=23)
	savefig(string("i_fr_", seq, ".png"), bbox_inches="tight", dpi=150)
	PyPlot.clf()
end


for (ind, val) in pairs(w8s)
	print(ind[1], ind[2], "\t", val, "\n")
	@show ind[1], ind[2]
end

totalinhHigh = zeros(4)
totalinhLow = zeros(4)
for ipop = 1:4
	emembers = filter(i->(i>0), popmembersHigh[ipop+12, :])
	totalinhHigh[ipop] = sum(weightsHigh[4001:4500, emembers])
	emembers = filter(i->(i>0), popmembersLow[ipop+4, :])
	totalinhLow[ipop] = sum(weightsLow[4001:4500, emembers])
end

plot(excfRates)
plot(inh2fRates)
plot(inh1fRates)
plot(restfRates)

plot(excfRates[26000:31000])
plot(inh2fRates[26000:31000])
plot(inh1fRates[26000:31000])
plot(restfRates[26000:31000])


pygui(true)
interval = 1000:5000
plot(succeSeqAct[1, interval])
plot(succeSeqAct[2, interval])
plot(succeSeqAct[3, interval])
plot(succeSeqAct[4, interval])
plot(succeSeqAct[5, interval])

# i2popmembers2 = findI2populations(w8si, Npop, popmembers2; iipop_len=25)
ipopmembers = findI2populations(new_weights, Npop, popmembers; iipop_len=25)
avg_weights = zeros(Npop, Npop)  # Contains the average strength from I-to-E
for iipop = 1:Npop
	imembers = ipopmembers[iipop,:]
	for ipop = 1:Npop
		members = filter(i->(i>0), popmembers[ipop, :])
		avg_weights[iipop, ipop] = sum(new_weights[imembers, members]) / count(i->(i>0), new_weights[imembers, members])
	end
end
x16 = reshape(avg_weights, (400, 1))
# x500_i = reshape(avg_weights, (400, 1))
x500_i
# dist = fit_mle(LogNormal, bins)
# (mu, sigma) = params(dist)
# sigma = std(x)
# y = ((1 ./ (sqrt(2 * pi) * sigma)) .* exp.(-0.5 .* (1 ./ sigma .* (bins .- mu)).^2)) # mu:mean, sigma:std (normal)
# y = (1 ./ bins .* sigma .* (sqrt(2*pi))) .* exp.(-(log.(bins .- mu) .^ 2) ./ 2*sigma^2) # (lognormal)
# plot(n)

# ipopmembers16 = findI2populations(weights16, Npop, popmembers16; iipop_len=25)

avg_weights = zeros(Npop, Npop)  # Contains the average strength from I-to-E
for iipop = 1:Npop
	imembers = ipopmembers[iipop,:]
	for ipop = 1:Npop
		members = filter(i->(i>0), popmembers[ipop, :])
		avg_weights[iipop, ipop] = sum(weightsOrig[imembers, members]) / count(i->(i>0), weightsOrig[imembers, members])
	end
end
xOrig = reshape(avg_weights, (400, 1))

#20, 28, 32, Orig, 40, 50


n, bins, patches = hist(x500, 30, density=true)
gaus_kernel = ImageFiltering.Kernel.gaussian((1.5,))
data = imfilter(n, gaus_kernel)
plot(bins, [data; 0])

pygui(true)


times = copy(timesHigh)
popmembers = copy(popmembersHigh)
ipopmembers = copy(ipopmembersHigh)
new_weights = copy(weightsHigh)

times = copy(timesLow)
popmembers = copy(popmembersLow)
ipopmembers = copy(ipopmembersLow)
new_weights = copy(weightsLow)

if true	# Single sequence
	pygui(false)
	Npop = size(popmembers, 1)
	Nmaxmembers = size(popmembers, 2)
	filter_min = 8050
	filter_max = 8300
	ipopsize = 25
	ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=ipopsize)
	seq = 5
	lineWidth = 1.2
	println("Creating plot...")
	figure(figsize=(3, 5))
	xlim(filter_min, filter_max)
	indexes = 1:seq_length
	indexes = indexes .+ ((seq-1)*seq_length)
	ylim(0, sum(popmembers[indexes, :] .> 0)+sum(ipopmembers[indexes, :] .> 0))#+500)
	ylabel(string("Sequence E"), fontsize="x-large")
	xlabel("Simulation Time (ms)", fontsize="x-large")
	tight_layout()
	ecolors = matplotlib.cm.Blues_r(1:round(Int,200/seq_length):200)
	icolors = matplotlib.cm.Reds_r(1:round(Int,200/seq_length):200)
	# Plot raster with the order of rows determined by population membership
	rowcount = 0
	ytickLoc = zeros(seq_length*2)
	for (ii, pp) in enumerate(indexes)
		print("\rpopulation ", pp)
		ecolor = ecolors[ii, :]
		for cc = 1:Nmaxmembers
			if popmembers[pp, cc] < 1
				break
			end
			global rowcount += 1
			ind = popmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			scatter(vals, y, s = .3, color=ecolor, marker="o", linewidths=lineWidth)
		end
		(mod(pp, seq_length) == 0) ? ind = seq_length : ind = mod(pp, seq_length)
		ytickLoc[ind] = rowcount
	end
	# Plot raster for interneurons (I2)
	for (ii, pp) in enumerate(indexes)
		icolor = icolors[ii, :]
		for cc = 1:ipopsize
			global rowcount += 1
			ind = ipopmembers[pp, cc]
			global vals = times[ind, 1:ns[ind]]
			global vals = filter(i->filter_max>i>filter_min, vals)
			global y = rowcount * ones(length(vals))
			scatter(vals, y, s = .3, color=icolor, marker="o", linewidths=lineWidth)
		end
		(mod(pp, seq_length) == 0) ? ind = seq_length + seq_length : ind = mod(pp, seq_length) + seq_length
		ytickLoc[ind] = rowcount
	end
	ax = gca()
	yticksL = zeros(length(ytickLoc))
	for (ii, tick) in enumerate(ytickLoc)
		if ii == 1
			yticksL[ii] = tick / 2
		else
			yticksL[ii] = ytickLoc[ii-1] + ((tick-ytickLoc[ii-1])/2)
		end
	end
	minor_ticks = matplotlib.ticker.FixedLocator(yticksL)

	ax.yaxis.set_major_locator(minor_ticks)
	ax.set_yticklabels([indexes; indexes], minor=false, fontsize=9, c="tab:gray")
	ax.yaxis.set_tick_params(which="minor", labelright=false, labelleft=false, right=false, left=false)
	ax.yaxis.set_tick_params(which="major", labelright=false, labelleft=true, right=false, left=false, direction="out", pad=1)
	for pp = 1:(seq_length*2)
		if mod(pp,seq_length) == 0
			ax.axhline([ytickLoc[pp]], linestyle="-", linewidth=.5, color="black")
			ax.axhline([ytickLoc[end]+(ipopsize * pp)], linestyle="-", linewidth=.3, color="black")
		else
			ax.axhline([ytickLoc[pp]], linestyle="-", linewidth=.15, color="black")
		end
	end
	ax.axhline([ytickLoc[seq_length]], linestyle="-", linewidth=.5, color="black")

	xticks = matplotlib.ticker.FixedLocator(filter_min:50:filter_max)
	ax.xaxis.set_major_locator(xticks)
	# ax.set_xticklabels(0:10:(filter_max-filter_min), minor=false, fontsize="x-large")
	ax.set_xticklabels(["0", "", "100", "", "200", ""], minor=false, fontsize="large")

	savefig(string("output_", T, ".png"), bbox_inches="tight", dpi=150)
	print("\rDone creating plot\n")
	PyPlot.clf()
end


seq_len = 4
Nseq = round(Int, Npop/seq_len)
NsuccAct = zeros(Npop)
NfailAct = zeros(Npop)
current = 0
expectation = 0

Nact = zeros(20)
for act in eAssemAct
	if act == 0
		continue
	else
		Nact[round(Int, act)] += 1
	end
end

for act in eAssemAct # Calculate the number of successful sequential dynamics
	if current == act || act == 0
		continue
	else    # First moment of assembly activation
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
		end
		(act in seq_len:seq_len:Npop) ? expectation = 0 : expectation = act + 1
	end
end

Npop = 20
avg_weights = zeros(Npop)
for pp = 1:Npop
	emembers = filter(i->(i>0), popmembers[pp, :])
	# for ip = 1:Npop
	# 	imembers = ipopmembers[ip, :]
	# 	avg_weights[ip, pp] = sum(new_weights[imembers, emembers]) / count(i->i>0,new_weights[imembers, emembers])
	# end
	avg_weights[pp] = sum(new_weights[4501:5000, emembers]) / count(i->i>0,new_weights[4501:5000, emembers])
end

avg_weightsE = zeros(Npop)
for pp = 1:Npop
	emembers = filter(i->(i>0), popmembers[pp, :])
	# for ip = 1:Npop
	# 	imembers = ipopmembers[ip, :]
	# 	avg_weights[ip, pp] = sum(new_weights[imembers, emembers]) / count(i->i>0,new_weights[imembers, emembers])
	# end
	# avg_weightsE[pp] = sum(new_weights[emembers, emembers]) / count(i->i>0,new_weights[emembers, emembers])
	avg_weightsE[pp] = sum(new_weights[1:4000, emembers]) / count(i->i>0,new_weights[1:4000, emembers])
end

x = zeros(20)
for pp = 1:20
	ip = pp + 1
	(ip == 21) && (ip=1)
	x[pp] = avg_weights[pp, ip]
end

i=2
for pp in diag(avg_weights)
	x[i] = pp
	i += 1
	(i==21) && (i=1)
end

eRates = copy(eRates20); popmembers = copy(popmembers20); ipopmembers = copy(ipopmembers20); new_weights = copy(w8s20)
eRates = copy(eRates50); popmembers = copy(popmembers50); ipopmembers = copy(ipopmembers50); new_weights = copy(w8s50)
eRates = copy(eRates100); popmembers = copy(popmembers100); ipopmembers = copy(ipopmembers100); new_weights = copy(w8s100)
eRates = copy(eRates200); popmembers = copy(popmembers200); ipopmembers = copy(ipopmembers200); new_weights = copy(w8s200)
eRates = copy(eRates500); popmembers = copy(popmembers500); ipopmembers = copy(ipopmembers500); new_weights = copy(w8s500)


excfRates20, restfRates20, inhfRates20, inh2fRates20 = getPopulationfRates(eRates20, Npop)
excfRates50, restfRates50, inhfRates50, inh2fRates50 = getPopulationfRates(eRates50, Npop)
excfRates100, restfRates100, inhfRates100, inh2fRates100 = getPopulationfRates(eRates100, Npop)
excfRates200, restfRates200, inhfRates200, inh2fRates200 = getPopulationfRates(eRates200, Npop)
excfRates500, restfRates500, inhfRates500, inh2fRates500 = getPopulationfRates(eRates500, Npop)


excfRates01, restfRates01, inhfRates01, inh2fRates01 = getPopulationfRates(eRates01, Npop)
excfRates05, restfRates05, inhfRates05, inh2fRates05 = getPopulationfRates(eRates05, Npop)
excfRates1, restfRates1, inhfRates1, inh2fRates1 = getPopulationfRates(eRates1, Npop)
excfRates10, restfRates10, inhfRates10, inh2fRates10 = getPopulationfRates(eRates10, Npop)
excfRates100, restfRates100, inhfRates100, inh2fRates100 = getPopulationfRates(eRates100, Npop)
excfRates500, restfRates500, inhfRates500, inh2fRates500 = getPopulationfRates(eRates500, Npop)


dataFULL = mean(avg_weights, dims=2)
dataNOIE = mean(avg_w8s, dims=2)
dataNOIE2 = mean(avg_w8s2, dims=2)

# dataFULL = mean(w8s[4501:5000, 1:4000], dims=2)
# dataNOIE = mean(new_weights[4501:5000, 1:4000], dims=2)

boxplot([dataFULL dataNOIE], vert=false, labels=["with eiSTDP", "without eiSTDP"])


avg_weights = zeros(Npop, Npop)
for pp = 1:Npop
	emembers = filter(i->(i>0), popmembers[pp, :])
	for ip = 1:Npop
		imembers = ipopmembers[ip, :]
		avg_weights[ip, pp] = sum(new_weights[imembers, emembers]) / count(i->i>0,new_weights[imembers, emembers])
	end
end
avg_weights2 = zeros(Npop, Npop)
for pp = 1:Npop
	emembers = filter(i->(i>0), popmembers2[pp, :])
	for ip = 1:Npop
		imembers = ipopmembers2[ip, :]
		avg_weights2[ip, pp] = sum(w8s[imembers, emembers]) / count(i->i>0,w8s[imembers, emembers])
	end
end
avg_weights3 = zeros(Npop, Npop)
for pp = 1:Npop
	emembers = filter(i->(i>0), popmembers3[pp, :])
	for ip = 1:Npop
		imembers = ipopmembers3[ip, :]
		avg_weights3[ip, pp] = sum(w8s2[imembers, emembers]) / count(i->i>0, w8s2[imembers, emembers])
	end
end
avg_weights4 = zeros(Npop, Npop)
for pp = 1:Npop
	emembers = filter(i->(i>0), popmembers4[pp, :])
	for ip = 1:Npop
		imembers = ipopmembers4[ip, :]
		avg_weights4[ip, pp] = sum(w8s3[imembers, emembers]) / count(i->i>0,w8s3[imembers, emembers])
	end
end

if true
	figure(figsize=(6, 6))
    xlim(0, 100)
    # ylim(0, 3.2)
	xlabel("Resting phase percentage (%)", fontsize="x-large")
    ylabel("Normalized Mean\nSynaptic Coupling (Pf)", fontsize="x-large")
	assembly = 1
	plot_dataEI = zeros(100)
	plot_dataIE = zeros(100)
	for perc = 1:100
		plot_dataEI[perc] = mean(avg_EI_tracker[perc])#avg_EI_tracker[perc][assembly, assembly+i]
		plot_dataIE[perc] = mean(avg_IE_tracker[perc])#avg_IE_tracker[perc][assembly, assembly+i]
	end
	plot(1:100, plot_dataEI ./ maximum(plot_dataEI), label="E→I₂", color="tab:blue", lw=2)
	plot(1:100, plot_dataIE ./ maximum(plot_dataIE), label="I₂→E", color="tab:red", lw=2)
	ax = gca()
	yticks = matplotlib.ticker.FixedLocator(0:0.2:1)
	ax.yaxis.set_major_locator(yticks)
	ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="large")
	ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], minor=false, fontsize="large")
	legend(fontsize="large") # loc="upper left", ncol=2, bbox_to_anchor=(0,1)
end

if true
	pygui(false)
	figure(figsize=(6, 4))
	ln_width = 2
    xlim(0,100)
    # ylim(55, 95)
    ylabel("Average I₂→E Assembly\nCoupling Strength (pF)", fontsize="xx-large")
    xlabel("Resting phase percentage (%)", fontsize="xx-large")
	pcolors = matplotlib.cm.Reds(56:round(Int, 200/Npop):256)
	plot_dataIE = zeros(100, Npop)
	for perc = 1:100
		for ipop = 1:Npop
			plot_dataIE[perc, ipop] = mean(avg_IE_tracker[perc][1:Npop, ipop])
		end
	end
	for ipop = 1:Npop
		plot(1:100, plot_dataIE[:, ipop], label="I₂→E", color=pcolors[ipop, :], lw=2)
	end
	ax = gca()
	ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
	grid(true, axis="y", ls="dotted")
	savefig(string("restAll_I2-E.png"), bbox_inches="tight")
    PyPlot.clf()
end
print()

if true
	pygui(false)
	figure(figsize=(6, 4))
	ln_width = 2
    xlim(0,100)
    # ylim(55, 95)
    ylabel("Average E→I₂ Assembly\nCoupling Strength (pF)", fontsize="xx-large")
    xlabel("Resting phase percentage (%)", fontsize="xx-large")
	pcolors = matplotlib.cm.Blues(56:round(Int, 200/Npop):256)
	plot_dataIE = zeros(100, Npop)
	for perc = 1:100
		for ipop = 1:Npop
			plot_dataEI[perc, ipop] = mean(avg_EI_tracker[perc][1:Npop, ipop])
		end
	end
	for ipop = 1:Npop
		plot(1:100, plot_dataEI[:, ipop], label="I₂→E", color=pcolors[ipop, :], lw=2)
	end
	ax = gca()
	ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
	grid(true, axis="y", ls="dotted")
	savefig(string("restAll_E-I2.png"), bbox_inches="tight")
    PyPlot.clf()
end
print()

if true
	pygui(false)
	figure(figsize=(6, 4))
	ln_width = 2
    xlim(0,100)
    # ylim(55, 95)
    ylabel("Average E & I₂ Assembly\nFiring Rates (Hz)", fontsize="xx-large")
    xlabel("Resting phase percentage (%)", fontsize="xx-large")
	ecolors = matplotlib.cm.Blues(56:10:256)
	icolors = matplotlib.cm.Reds(56:10:256)
	plot_datafr = fRates
	# plot_datafr = zeros(100, 44)
	# for perc = 1:100
	# 	for ipop = 1:44
	# 		plot_datafr[perc, ipop] = mean(fr_tracker[perc][ipop])
	# 	end
	# end
	for ipop = 1:20#21:40
		plot(2:100, plot_datafr[ipop, 2:end], label="I₂→E", color=ecolors[ipop, :], lw=2)
	end
	for ipop = 21:40#21:40
		plot(2:100, plot_datafr[ipop, 2:end], label="I₂→E", color=icolors[ipop-20, :], lw=2)
	end
	# plot(2:100, plot_datafr[41, 2:end], label="Excitatory", color="tab:blue", lw=2)
	# plot(2:100, plot_datafr[42, 2:end], label="Inhibitory", color="tab:orange", lw=2)
	# plot(2:100, plot_datafr[43, 2:end], label="I₁", color="tab:brown", lw=2)
	# plot(2:100, plot_datafr[44, 2:end], label="I₂", color="tab:red", lw=2)
	ax = gca()
	ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="x-large")
    ax.xaxis.set_tick_params(which="major", labeltop=false, labelbottom=true, top=false, bottom=false, labelsize="x-large")
	grid(true, axis="y", ls="dotted")
	savefig(string("restAll_fr.png"), bbox_inches="tight")
    PyPlot.clf()
end
print()

# if true	# I2-to-E recurrent coupling comparison
# 	figure(figsize=(20, 6))
#     # xlim(0, 100)
#     # ylim(0, 3.2)
# 	xlabel("Resting phase percentage (%)", fontsize="x-large")
#     ylabel("Average E/I₂ Assembly\nSynaptic Coupling (Pf)", fontsize="x-large")
# 	assembly = 1
# 	i=15
# 	plot_dataEI = zeros(100, 20)
# 	plot_dataIE = zeros(100, 19)
# 	for perc = 1:100
# 		# plot_dataEI[perc, :] = avg_EI_tracker[perc][diagind(avg_EI_tracker[perc])]#avg_EI_tracker[perc][assembly, assembly+i]
# 		# plot_dataIE[perc, :] = avg_IE_tracker[perc][diagind(avg_IE_tracker[perc])]#avg_IE_tracker[perc][assembly, assembly+i]
# 		plot_dataIE[perc, :] = avg_IE_tracker[perc][1:Npop, assembly]
# 	end
# 	x = 1:100
# 	# bar(x, plot_dataIE[100, :], color="tab:blue", width=0.5, label="After resting phase")
# 	# bar((x.+0.2), plot_dataIE[1, :], color="tab:cyan", width=0.5, label="Before resting phase")
#     # [box.set_facecolor("tab:blue") for box in plt["boxes"]]
# 	bar(x, plot_dataIE[100, :], color="tab:blue", width=0.5, label="After resting phase")
# 	bar((x.+0.2), plot_dataIE[1, :], color="tab:cyan", width=0.5, label="Before resting phase")
# 	ax = gca()
# 	x_ticks = matplotlib.ticker.FixedLocator(1:Npop)
# 	ax.xaxis.set_major_locator(x_ticks)
# 	ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="large")
# 	ax.yaxis.set_tick_params(which="major", labelleft=true, labelright=false, left=true, right=false, labelsize="large")
# 	legend(fontsize="large") # loc="upper left", ncol=2, bbox_to_anchor=(0,1)
# end

if true	# I2-to-E average develop during rest
	# figure(figsize=(5, 5))
	_, axs = subplots(2, 1, sharex=true)
	xlim(0, 100)
	xlabel("Resting phase percentage (%)", fontsize="x-large")
    ylabel("Average E/I₂ Assembly\nSynaptic Coupling (Pf)", fontsize="x-large")
	plot_dataEI = zeros(100)
	plot_dataIE = zeros(100)
	plot_dataEI_diag = zeros(100)
	plot_dataIE_diag = zeros(100)
	plot_dataEI_noDiag = zeros(100)
	plot_dataIE_noDiag = zeros(100)
	perc = 1
	for iter = 1:20:1981
		plot_dataEI[perc] = mean(avgEI[:, iter:(iter+19)])
		plot_dataIE[perc] = mean(avgIE[:, iter:(iter+19)])
		plot_dataEI_diag[perc] = mean(diag(avgEI[:, iter:(iter+19)]))
		plot_dataIE_diag[perc] = mean(diag(avgIE[:, iter:(iter+19)]))
		plot_dataEI_noDiag[perc] = (sum(avgEI[:, iter:(iter+19)]) - sum(diag(avgEI[:, iter:(iter+19)]))) / 380
		plot_dataIE_noDiag[perc] = (sum(avgIE[:, iter:(iter+19)]) - sum(diag(avgIE[:, iter:(iter+19)]))) / 380
		perc += 1
	end
	plot_dataEI ./= maximum(plot_dataEI)
	# plot_dataIE ./= maximum(plot_dataIE)
	plot_dataEI_diag ./= maximum(plot_dataEI_diag)
	# plot_dataIE_diag ./= maximum(plot_dataIE_diag)
	plot_dataEI_noDiag ./= maximum(plot_dataEI_noDiag)
	x = 1:100
	axs[1].plot(x, plot_dataEI, label="Avg E-to-I", color="tab:blue")
	axs[2].plot(x, plot_dataIE, label="Avg I-to-E", color="tab:red")
	axs[1].plot(x, plot_dataEI_diag, label="Avg rec. E-to-I", color="blue")
	axs[2].plot(x, plot_dataIE_diag, label="Avg rec. I-to-E", color="red")
	axs[1].plot(x, plot_dataEI_noDiag, label="Avg no rec. E-to-I", color="darkblue")
	axs[2].plot(x, plot_dataIE_noDiag, label="Avg no rec. I-to-E", color="brown")
	axs[1].legend(fontsize="large") # loc="upper left", ncol=2, bbox_to_anchor=(0,1)
	axs[2].legend(fontsize="large")
	savefig("averageWeights_rest.png", bbox_inches="tight", dpi=150)
	PyPlot.clf()
end

if true	# I2-to-E recurrent average develop during rest
	figure(figsize=(5, 5))
	xlabel("Resting phase percentage (%)", fontsize="x-large")
    ylabel("Average E/I₂ Assembly\nSynaptic Coupling (Pf)", fontsize="x-large")
	plot_dataEI = zeros(100)
	plot_dataIE = zeros(100)
	perc = 1
	for iter = 1:20:1981
		plot_dataEI[perc] = mean(diag(avgEI[:, iter:(iter+19)]))
		plot_dataIE[perc] = mean(diag(avgIE[:, iter:(iter+19)]))
		perc += 1
	end
	plot_dataEI ./= maximum(plot_dataEI)
	plot_dataIE ./= maximum(plot_dataIE)
	x = 1:100
	plot(x, plot_dataEI, label="Avg E-to-I")
	plot(x, plot_dataIE, label="Avg I-to-E")
	ax = gca()
	legend(fontsize="large") # loc="upper left", ncol=2, bbox_to_anchor=(0,1)
	savefig("averageWeights_rest.png", bbox_inches="tight", dpi=150)
	PyPlot.clf()
end

if true	# I2-to-E recurrent weights
	# figure(figsize=(5, 5))
	_, axs = subplots(2, 1, sharex=false)
	seq = 1
	nbins = 70
	flt_window=2.5
	ln_width = 2.
	alpha = .3
	hcolors1 = matplotlib.cm.Blues_r(100:20:200)
	pcolors1 = copy(hcolors1)
    hcolors1[:, 4] .= alpha
	hcolors2 = matplotlib.cm.Reds_r(100:20:200)
	pcolors2 = copy(hcolors2)
    hcolors2[:, 4] .= alpha
	xlabel("Resting phase percentage (%)", fontsize="x-large")
    ylabel("Average E/I₂ Assembly\nSynaptic Coupling (pF)", fontsize="x-large")
	plot_dataEI = zeros(5, 20)
	plot_dataIE = zeros(5, 20)
	perc = 1
	for iter = 1:500:1981
		plot_dataEI[perc, :, :] = diag(avgEI[:, iter:(iter+19)])
		plot_dataIE[perc, :, :] = diag(avgIE[:, iter:(iter+19)])
		perc += 1
	end
	plot_dataEI[perc, :, :] = diag(avgEI[:, 1981:(1981+19)])
	plot_dataIE[perc, :, :] = diag(avgIE[:, 1981:(1981+19)])

	plt1 = axs[1].boxplot([plot_dataEI[1, :], plot_dataEI[2, :], plot_dataEI[3, :], plot_dataEI[4, :], plot_dataEI[5, :]], labels=["0%", "25%", "50%", "75%", "100%"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2, "color"=>"tab:orange"), sym="")
	[box.set_facecolor("tab:blue") for box in plt1["boxes"]]
	plt2 = axs[2].boxplot([plot_dataIE[1, :], plot_dataIE[2, :], plot_dataIE[3, :], plot_dataIE[4, :], plot_dataIE[5, :]], labels=["0%", "25%", "50%", "75%", "100%"], patch_artist=true, medianprops=Dict("linestyle"=>"-", "linewidth"=>2, "color"=>"tab:cyan"), sym="")
	[box.set_facecolor("tab:red") for box in plt2["boxes"]]
	# axs[1].legend(fontsize="large") # loc="upper left", ncol=2, bbox_to_anchor=(0,1)
	# axs[2].legend(fontsize="large")
	savefig("boxWeights_rest.png", bbox_inches="tight", dpi=150)
	PyPlot.clf()
end
restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])


fid = h5open("C:/Users/disma/Desktop/MSc Thesis/Code/litwin-kumar_doiron_formation_2014/data_test_NoOver_2i.h5", "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
close(fid)

avgE = zeros(20)
avgE_new = zeros(20)
avgIE = zeros(20)
avgIE_new = zeros(20)
avgEI = zeros(20)
avgEI_new = zeros(20)
for ipop = 1:Npop
	emembers = filter(i->(i>0), popmembers[ipop, :])
	imembers = ipopmembers[ipop, :]
	avgE[ipop] = sum(weights[emembers, emembers]) / count(i->(i>0), weights[emembers, emembers])
	avgE_new[ipop] = sum(new_weights[emembers, emembers]) / count(i->(i>0), new_weights[emembers, emembers])
	avgIE[ipop] = sum(weights[imembers, emembers]) / count(i->(i>0), weights[imembers, emembers])
	avgIE_new[ipop] = sum(new_weights[imembers, emembers]) / count(i->(i>0), new_weights[imembers, emembers])
	avgEI[ipop] = sum(weights[emembers, imembers]) / count(i->(i>0), weights[emembers, imembers])
	avgEI_new[ipop] = sum(new_weights[emembers, imembers]) / count(i->(i>0), new_weights[emembers, imembers])
end

# --- LOAD AND SAVE DATA ---
#____________________________________________________________________________

# NOTE:: DO NOT DELETE THIS TO REMEMBER HOW EVERYTHING IS SAVED!!!
fid = h5open("saved_Qvals.h5","w")
g = create_group(fid, "data")
g["xValsTau"] = [x20 x50 xOrig x200 x500]
g["xValsiLamda"] = [x01 x05 xOrig x10 x100 x500_i]
g["qValsTau"] = [q20 q50 qOrig q200 q500_t]
g["qValsiLamda"] = [q01 q05 qOrig q10 q100 q500_i]
close(fid)


# NOTE:: DO NOT DELETE THIS TO REMEMBER HOW EVERYTHING IS SAVED!!!
fid = h5open("saved_IIvals.h5","w")
g = create_group(fid, "data")
g["eRates16"] = eRates16
g["eRates20"] = eRates20
g["eRates28"] = eRates28
g["eRates32"] = eRates32
g["eRates40"] = eRates40
g["eRates50"] = eRates50
g["popmembers"] = [popmembers16 popmembers20 popmembers28 popmembers32 popmembers40 popmembers50]
g["ipopmembers"] = [ipopmembers16 ipopmembers20 ipopmembers28 ipopmembers32 ipopmembers40 ipopmembers50]
g["weights"] = [w8s16 w8s20 w8s28 w8s32 w8s40 w8s50]
close(fid)


# TO LOAD THE DATA:
fid = h5open("saved_IIvals.h5", "r")
eRates16 = read(fid["data"]["eRates16"])
eRates20 = read(fid["data"]["eRates20"])
eRates28 = read(fid["data"]["eRates28"])
eRates32 = read(fid["data"]["eRates32"])
eRates40 = read(fid["data"]["eRates40"])
eRates50 = read(fid["data"]["eRates50"])
popmembersL = read(fid["data"]["popmembers"])
ipopmembers = read(fid["data"]["ipopmembers"])
weights = read(fid["data"]["weights"])
close(fid)


ipopmembers16 = ipopmembers[:, 1:25]
ipopmembers20 = ipopmembers[:, 26:50]
ipopmembers28 = ipopmembers[:, 51:75]
ipopmembers32 = ipopmembers[:, 76:100]
ipopmembers40 = ipopmembers[:, 101:125]
ipopmembers50 = ipopmembers[:, 126:150]

popmembers16 = popmembersL[:, 1:300]
popmembers20 = popmembersL[:, 301:600]
popmembers28 = popmembersL[:, 601:900]
popmembers32 = popmembersL[:, 901:1200]
popmembers40 = popmembersL[:, 1201:1500]
popmembers50 = popmembersL[:, 1501:1800]

weights16 = weights[:, 1:5000]
weights20 = weights[:, 5001:10000]
weights28 = weights[:, 10001:15000]
weights32 = weights[:, 15001:20000]
weights40 = weights[:, 20001:25000]
weights50 = weights[:, 25001:30000]




fid = h5open("saved_IIvals_extended.h5","w")
g = create_group(fid, "data")
g["eRatesOrig"] = eRatesOrig
g["popmembersOrig"] = popmembersOrig
g["ipopmembersOrig"] = ipopmembersOrig
g["weightsOrig"] = weightsOrig
g["avg_inter"] = [avg_inter16 avg_inter20 avg_inter28 avg_interOrig avg_inter32 avg_inter40 avg_inter50]
g["qVals"] = [q16 q20 q28 qOrig q32 q40 q50]
g["total_inter"] = [total_inter16 total_inter20 total_inter28 total_interOrig total_inter32 total_inter40 total_inter50]
g["xVals"] = [x16 x20 x28 xOrig x32 x40 x50]
g["yVals"] = [y16 y20 y28 yOrig y32 y40 y50] # For plotting population firing rates
close(fid)

# Last from rest
saved_dataEI = avg_EI_tracker[1]
saved_dataIE = avg_IE_tracker[1]
saved_datafr = fr_tracker[1]
for perc = 2:100
	saved_dataEI = [saved_dataEI avg_EI_tracker[perc]]
	saved_dataIE = [saved_dataIE avg_IE_tracker[perc]]
	saved_datafr = [saved_datafr fr_tracker[perc]]
end
fid = h5open("saved_restVals.h5","w")
g = create_group(fid, "data")
g["popmembers"] = popmembers
g["ipopmembers"] = ipopmembers
g["firing_rates"] = saved_datafr
g["avg_EI"] = saved_dataEI
g["avg_IE"] = saved_dataIE
close(fid)


fid = h5open("saved_restVals.h5","r")
popmembers = read(fid["data"]["popmembers"])
ipopmembers = read(fid["data"]["ipopmembers"])
fRates = read(fid["data"]["firing_rates"])	# [eAssemfRates; iAssemfRates; excfRates; inhfRates; inh1fRates; inh2fRates]
avgEI = read(fid["data"]["avg_EI"])
avgIE = read(fid["data"]["avg_IE"])
close(fid)

















# REALLY DEBUGGING:

fig = figure(figsize=(25, 5), constrained_layout=true)
# fig = figure(tight_layout=true)
plotRange = 1:1000
linewidth = .5

gs = matplotlib.gridspec.GridSpec(5, 25, figure=fig)

spikes, epsc, ipsc, epsp, ipsp = runSingleSynapse() 

ax1 = subplot(gs.new_subplotspec((1, 1), colspan=12))
ax1.plot(plotRange, spikes[plotRange], lw=linewidth)
ax1.plot(plotRange, epsc[plotRange] .- 7, lw=linewidth)
ax1.set_title("Excitatory synapse")
ax1.tick_params(left=false, right=false, labelleft=false, labelright=false, labelbottom=false)

ax2 = subplot(gs.new_subplotspec((2, 1), colspan=12, rowspan=4))
ax2.plot(plotRange, epsp[plotRange], lw=linewidth)
ax2.set(ylabel="Post-synaptic potential (mV)")
ax2.set(xlabel="Time (ms)")
ax2.tick_params()

ax3 = subplot(gs.new_subplotspec((1, 13), colspan=12))
ax3.plot(plotRange, spikes[plotRange], lw=linewidth)
ax3.plot(plotRange, ipsc[plotRange] .- 7, lw=linewidth)
ax3.set_title("Inhibitory synapse")
ax3.tick_params(left=false, right=false, labelleft=false, labelright=false, labelbottom=false)

ax4 = subplot(gs.new_subplotspec((2, 13), colspan=12, rowspan=4))
ax4.plot(plotRange, ipsp[plotRange], lw=linewidth)
ax4.set(xlabel="Time (ms)")
ax4.tick_params()







eneurons = filter(i->(i>0), unique(popmembers))
ineurons = unique(ipopmembers)
if true
	# efrates = 0
	# ifrates = 0
	# for pp = 1:Npop
	# 	ifrates += count(i->(i>0), times[ipopmembers[pp, :]])
	# 	efrates += count(i->(i>0), times[filter(i->(i>0), popmembers[pp, :])])
	# end
	efrates = count(i->(i>0), times[eneurons, :])
	ifrates = count(i->(i>0), times[ineurons, :])
	efrates = efrates / length(eneurons)
	ifrates = ifrates / length(ineurons)

	efrates = 1000 * efrates / T
	ifrates = 1000 * ifrates / T
end
