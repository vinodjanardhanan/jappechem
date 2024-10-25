using GlobalSensitivity, QuasiMonteCarlo, Statistics, Plots, CSV, DataFrames, Printf
export eis_gsa, eisfit

#function to read the elechemical model parameters from CSV file and returns the parameters as Array
function read_ecm_model_parameters()

	df = CSV.read("ecm.csv", DataFrame)
	params = Array{Float64,1}()
	for i in eachrow(df)
		for j in i
			push!(params,j)
		end
	end
	params
end

# function for doing the global sensitivity analysis
function eis_gsa()
	#create output file for Sobol analysis
	sobol_S1 = open("sobol_S1.csv", "w")
	sobol_ST = open("sobol_ST.csv", "w")
	@printf(sobol_S1,"%s","R0,R1,R2,R3,Q1,Q2,Q3,a1,a2,a3,WR,WT,WP\n")
	@printf(sobol_ST,"%s","R0,R1,R2,R3,Q1,Q2,Q3,a1,a2,a3,WR,WT,WP\n")
	samples = 15000
	p = read_ecm_model_parameters()
	lb = (1-0.1) .* p
	ub = (1+0.1) .* p
	sampler = SobolSample()
	A, B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)

	freq = freq_settings()
	
	function eis(params)
		ecm_fit(freq, params; sens = true)
	end
	
	sobol_indices = gsa(eis, Sobol(), A, B)
	S1 = string.(sobol_indices.S1[1,:])
	ST = string.(sobol_indices.ST[1,:])
	@printf(sobol_S1, "%s", join(S1,","))
	@printf(sobol_ST, "%s", join(ST,","))
	close(sobol_S1)
	close(sobol_ST)
	
end


#function for reading the frequency related information if EIS experiment from CSV file
function freq_settings()
	df = CSV.read("freqdata.csv", DataFrame)
	ω_min = df.w_min[1]
	ω_max = df.w_max[1]
	n_dec = ceil(log10(ω_max/ω_min))
	step_size = 10^(1/df.n_step[1])
	ω =  ω_min
	ωs = Array{Float64,1}()
	for i in 1:n_dec	
		for k in 1:df.n_step[1]
			push!(ωs, ω)
			ω *= step_size
		end
	end
	ωs
end

#function for performing the EIS fit using the ECM model parameters
function ecm_fit(ωs, params; sens = false)
	
	#open a file for writing output, this file will be used only if sens = false
	ecmout = open("modeloutput.csv","w")
	z_real = Array{Float64,1}()
	z_imag = Array{Float64,1}()
	zmag = Array{Float64, 1}()
	R0, R1, R2, R3 = params[1:4]
	Q1, Q2, Q3 = params[5:7]
	a1, a2, a3 = params[8:10]
	WR, WT, WP = params[11:end]

	for ω in ωs
		z_rcp1 = R1/(1+(ω*im)^a1 * Q1 * R1)
		z_rcp2 = R2/(1+(ω*im)^a2 * Q2 * R2)
		z_rcp3 = R3/(1+(ω*im)^a3 * Q3 * R3)
		z_w = WR*tanh((WT*ω*im)^WP)/(WT*ω*im)^WP
		z_total = R0 + z_rcp1 + z_rcp2 + z_rcp3 + z_w
		zr = real(z_total)
		zi = imag(z_total)
		push!(zmag, abs(z_total))
		push!(z_real, zr)
		push!(z_imag, zi)
		data = string(zr)*","* string(zi)
		@printf(ecmout, "%s", data)
		@printf(ecmout, "\n")
	end
	close(ecmout)
	if sens		
		zmag
	end
	
end

function ecm_fit()
	params = read_ecm_model_parameters()
	freq = freq_settings()
	ecm_fit(freq, params)
end
