function do_subplot(x, y, z)
############################    
    ax = subplot(x, y, z, polar="true") 
    ax.grid(:true, fontsize=:10)
    ax.set_rlabel_position(-90)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    PyPlot.title(string(z))

end    # do_subplot()
    

# smooth the spectra into bands centered on 0.05Hz spacing (i.e. 0:0.005:0.64)
function smooth_spectra(Pden_in, sample_frequency)
##################################################

    nyquist = sample_frequency/2

    freq_in = range(0, stop=nyquist, length=length(Pden_in))

    freq_out = [0.0]
    Pden_smoothed = [mean(Pden_in[1:8])]

    i = 9
    while i <= length(Pden_in)

        push!(freq_out,freq_in[i+8])

        if i < length(Pden_in)-16

            push!(Pden_smoothed, mean(Pden_in[i:i+16]))

        end

        i+=16

    end

    push!(Pden_smoothed, mean(Pden_in[end-8:end]))
            
    return (freq_out, Pden_smoothed)
        
end


function get_Fourier_coefficients(heave, north, west)
#####################################################  
    
    # Get the cross periodograms
    cps_heave_heave = mt_cross_power_spectra([heave heave]', fs=sample_frequency);
    cps_north_north = mt_cross_power_spectra([north north]', fs=sample_frequency);
    cps_west_west = mt_cross_power_spectra([west west]', fs=sample_frequency);

    cps_north_heave = mt_cross_power_spectra([north heave]', fs=sample_frequency);
    cps_west_heave = mt_cross_power_spectra([west heave]', fs=sample_frequency);
    cps_north_west = mt_cross_power_spectra([north west]', fs=sample_frequency);

##    fhh = cps_heave_heave.freq
    fhh, Chh = smooth_spectra(real.(cps_heave_heave.power[1,1,:]), sample_frequency)

    #fnn = cps_north_north.freq
    fhh, Cnn = smooth_spectra(real.(cps_north_north.power[1,1,:]), sample_frequency)

    #fww = cps_west_west.freq
    fhh, Cww = smooth_spectra(real.(cps_west_west.power[1,1,:]), sample_frequency)

    #fnw = cps_north_west.freq
    fhh, Cnw = smooth_spectra(real.(cps_north_west.power[1,2,:]), sample_frequency)

    #fnh = cps_north_heave.freq
    fhh, Qnh = smooth_spectra(imag.(cps_north_heave.power[1,2,:]), sample_frequency)

    #fwh = cps_west_heave.freq
    fhh, Qwh = smooth_spectra(imag.(cps_west_heave.power[1,2,:]), sample_frequency)

    a1 = Qnh ./ ((Cnn .+ Cww) .* Chh) .^ 0.5
    b1 = -Qwh ./ ((Cnn .+ Cww) .* Chh) .^ 0.5

    a2 = (Cnn .- Cww) ./ (Cnn .+ Cww)
    b2 = -2 .* Cnw ./ (Cnn .+ Cww)
    
    return(fhh, Chh, a1, b1, a2, b2)
    
end    # get_Fourier_coefficients()


function get_displacements(arry)
#####################################
    
    displacements = []

    if length(arry[1]) == 3
    
        for i ∈ arry
            append!(displacements,parse(Int, SubString.(i, 1, 1), base=16)*16^2 + parse(Int, SubString.(i, 2, 2), base=16)*16^1 + parse(Int, SubString.(i, 3, 3), base=16)*16^0)
        end
        
    else
        
        for i ∈ arry
            append!(displacements,parse(Int, SubString.(i, 1, 1), base=16)*16^1 + parse(Int, SubString.(i, 2, 2), base=16)*16^0)
        end
        
    end

    displacements[findall(>=(2048), displacements)] = 2048 .- displacements[findall(>=(2048), displacements)];
    
    return(displacements./100)
    
end     # get_displacements()


function get_HNW(infil)
#####################################
        
    global df = DataFrame(CSV.File(infil,header=0, delim=",", types=String));

    # Calculate sequence numbers
    arry = SubString.(df.Column1, 3, 4)

    global sequence = []

    for i ∈ arry
        append!(sequence,parse(Int, SubString.(i, 1, 1), base=16)*16^1 + parse(Int, SubString.(i, 2, 2), base=16)*16^0)
    end

    arry = SubString.(df.Column3, 1, 3);
    heave = get_displacements(arry);

    # Calculate north WSEs
    arry = SubString.(df.Column3, 4, ) .* SubString.(df.Column4, 1, 2)
    north = get_displacements(arry);

    # Calculate north WSEs
    arry = SubString.(df.Column4, 3, 4) .* SubString.(df.Column5, 1, 1)
    west = get_displacements(arry);

    return(heave, north, west)

    end    # get_HNW()


function plot_polar(fig, displacement_df, row, col, total, spec_max)
####################################################################

# Called by: Do_POLAR_PLOTS

    Chh = displacement_df.Chh[total]
    a1 = displacement_df.a1[total]
    b1 = displacement_df.b1[total] 
    a2 = displacement_df.a2[total] 
    b2 = displacement_df.b2[total]
    time_string  = displacement_df.Time_string[total]

    aa = length(Chh) # Number of spectral points

    r = 1:6:aa
    global ρ = r ./ (aa/nyquist) 

    global θ = 0:pi/180:2pi        

    # populate a matrix of spectral surface values
    global mat = [Chh[r1] * (a1[r1]*cos.(θ) + b1[r1]*sin.(θ) + a2[r1]*cos.(2θ) + b2[r1]*sin.(2θ)) for r1 ∈ r]

    mat = hvcat(size(mat,1),mat...)

    # set any values less than zero to zero
    global mat[mat .< 0] .= 0

    time_string = string(total)*" - "*displacement_df.Time_string[total]

    cmap = Reverse(:ocean)
    levels = round(spec_max/100, digits=2):round(spec_max/20, digits=2):round(spec_max, digits=2)

    ax = CairoMakie.PolarAxis(fig[row, col],
    thetaticklabelsize = 15,  
    rlimits=(0,0.4), rticklabelsize=15, rticks=0:0.2:0.4, rgridwidth=0.5, rtickangle=180, rminorgridvisible=true, rminorgridstyle=:dot,
    theta_0=-pi/2, thetagridwidth=0.5, thetaminorgridvisible=true, thetaminorgridstyle=:dot, thetaminorticks=IntervalsBetween(3), 
    direction=-1, width=330, height=310, title=time_string, titlesize=18,
    )

    CairoMakie.contourf!(ax, θ, ρ, Float64.(mat), colormap=cmap, levels=levels) # 0.02:0.05:1,)
    CairoMakie.contour!(ax, θ, ρ, Float64.(mat), colormap=cmap, levels=levels)
    
    return(fig)
    
    end    # plot_polar()


function plot_spectra_and_direction(fig, displacement_df, row, col, total, spec_max)
####################################################################

# Called by Do_Spectral_Plots

    global fhh = displacement_df.fhh[total]
    Chh = displacement_df.Chh[total]
    a1 = displacement_df.a1[total]
    b1 = displacement_df.b1[total] 
    a2 = displacement_df.a2[total] 
    b2 = displacement_df.b2[total]
    f2 = displacement_df.f2[total]
    Pden2 = displacement_df.Pden2[total]
    time_string  = string(total)*" - "*displacement_df.Time_string[total]

    # Calculate the Centred Fourier Coefficients
    θ₀ = atan.(b1,a1)
    m1 = (a1.^2 .+ b1.^2).^0.5
##    m2 = a2.*cos.(2*θ₀) .+ b2.*sin.(2*θ₀)
##    n2 = -a2.*sin.(2*θ₀) .+ b2.*cos.(2*θ₀)

    # Calculate the spread
    σc = (2 .* (1 .- m1)).^0.5;

    direction = mod.(rad2deg.(atan.(b1,a1)),360)
    
    ###################################################################################

    function nearest_neighbor(target, points)
        
        # use a nearest neighbour approach to locate the closest point
        distances = [norm(target - point) for point in points]
        nearest_index = argmin(distances)
        return points[nearest_index]
    end

    new_dir = []

    current = direction[1]

    for next ∈ direction[2:end]

        push!(new_dir,current)
        
        # create three points to cover full range from 0-360⁰
        points = [next - 360, next, next + 360]

        # now find the point nearest to current point
        nearest = nearest_neighbor(current, points)
        current = nearest

    end

    push!(new_dir,current)
    new_dir = Float32.(new_dir)
    
    ###################################################################################

    ax1 = fig[row,col] = Axis(fig, limits=(0, 0.64, 0, 360),
        xlabel="Frequency (Hertz)", ylabel = "Direction (⁰)", yreversed=true, yticks=0:45:360, width=330, height=310, title=time_string, titlesize=18)
    ax1.xlabelsize=20; ax1.ylabelsize=20
        
    ax2 = fig[row,col] = Axis(fig,  limits=(0, 0.64, 0, spec_max),
        ylabel="Spectral Density (m²/Hz.)", yaxisposition=:right, yticklabelalign=(:left, :center), ylabelrotation=-pi/2)
    ax2.xlabelsize=20; ax2.ylabelsize=20
        
    hidedecorations!(ax2, ticks=false, ticklabels=false, label=false)

    CairoMakie.lines!(ax1, fhh, new_dir, linewidth=0.75, color=:grey)
    CairoMakie.lines!(ax1, fhh, new_dir.+360, linewidth=0.75, color=:grey)
    CairoMakie.lines!(ax1, fhh, new_dir.-360, linewidth=0.75, color=:grey)

    CairoMakie.band!(ax1, fhh, new_dir .+ rad2deg.(σc), new_dir .- rad2deg.(σc), fillrange = new_dir .- rad2deg.(σc), color=(:grey, 0.125), )
    CairoMakie.band!(ax1, fhh, new_dir .+360 .+ rad2deg.(σc), new_dir .+360 .- rad2deg.(σc), fillrange = new_dir .- rad2deg.(σc), color=(:grey, 0.125), )
    CairoMakie.band!(ax1, fhh, new_dir .-360 .+ rad2deg.(σc), new_dir .-360 .- rad2deg.(σc), fillrange = new_dir .- rad2deg.(σc), color=(:grey, 0.125), )

    CairoMakie.lines!(ax2, f2, Pden2, color=:red, linewidth=:2, fillrange=0, fillalpha=0.75, fillcolor=:red, z_order=:2)
            
    return(fig)
    
    end    # plot_spectra_and_direction()


# function to calc spectral moments
function calc_moments(t,y)
#########################
    
"""
    Calls: Nil
    Called by: calc_PM_JONSWAP()
"""
    
#==    
    moments = []    
    
    [push!(moments,sum([t[i]^n .* y[i] * 0.00125 for i ∈ 1:length(y)])) for n ∈ [-1,0,1,2,4]]

    m₋₁ = moments[1]; m₀ = moments[2]; m₁ = moments[3]; m₂ = moments[4]; m₄ = moments[5]
==#        
    ax1 = (last(t) - first(t)) / (length(t)-1)

    # calc spectral moments m0, m1, m2, m3, and m4
    s_01 = 0; s00 = 0; s01 = 0; s02 = 0; s03 = 0; s04 = 0;
    m₋₁ = 0; m₀ = 0; m₁ = 0; m₂ = 0; m₃ = 0; m₄ = 0

    for ii in 1:length(t)

        s_01 += t[ii]^-1 * y[ii]
        s00 += t[ii]^0 * y[ii]
        s01 += t[ii]^1 * y[ii]
        s02 += t[ii]^2 * y[ii]
        s03 += t[ii]^3 * y[ii]
        s04 += t[ii]^4 * y[ii]

    end

    m₋₁ = 0.5*ax1*(first(t)^-1*first(y) + 2*s_01 + last(t)^0*last(y))
    m₀ = 0.5*ax1*(first(t)^0*first(y) + 2*s00 + last(t)^0*last(y))
    m₁ = 0.5*ax1*(first(t)^1*first(y) + 2*s01 + last(t)^1*last(y))
    m₂ = 0.5*ax1*(first(t)^2*first(y) + 2*s02 + last(t)^2*last(y))
    m₃ = 0.5*ax1*(first(t)^3*first(y) + 2*s03 + last(t)^3*last(y))
    m₄ = 0.5*ax1*(first(t)^4*first(y) + 2*s04 + last(t)^4*last(y))        

    return(m₋₁, m₀, m₁, m₂, m₄)

    end    # calc_moments()
        

# function to calculate PM and JONSWAP spectra
function calc_PM_JONSWAP(frequency, spectra)
############################################

    α = 0.0081

    fₚ = frequency[argmax(spectra)]
    Tₚ = 1/fₚ

    m₋₁, m₀, m₁, m₂, m₄ = calc_moments(frequency,spectra)
    Hₘ₀ = 4 * √(m₀)
    A = 0.039370190176622376
    αⱼ = α #A * Hₘ₀^2 * Tₚ^-4 # Normalization factor
    σ = [f <= fₚ ? 0.07 : 0.09 for f ∈ frequency]
    r = [exp(-1 * ((f - fₚ)^2) / ((2 * σ[i]^2) * fₚ^2)) for (i, f) ∈ enumerate(frequency)]

    PM = [α * g^2 * (2π)^-4 * ff^-5 * exp(-(5/4) * (fₚ/ff)^4) for ff ∈ frequency]    # from Tucker and Pitt (2001) 5.5-3 p.100
    PM = [x == 0.0 ? NaN : x for x in PM]

    JONSWAP = [αⱼ * f^-5 * exp(-1.25 * ((f/fₚ)^-4)) * γ^r[i] for (i, f) ∈ enumerate(frequency)]

    return(PM, JONSWAP)

end    # calc_PM_JONSWAP()


# function to eliminate minor peaks where the trough between adjacent peaks is greater than 50% of the smaller peak    
function eliminate_minor_peaks(x_peaks, y_peaks, x_troughs, y_troughs)
######################################################################
 
"""
    Calls: Nil
    Called by: find_peaks_and_valleys()
"""
    
    # Initialize new arrays to store significant peaks and troughs
    x_peaks_significant = Float64[]
    y_peaks_significant = Float64[]
    x_troughs_significant = Float64[]
    y_troughs_significant = Float64[]

    # Loop over troughs (the number of veal troughs is 1 less than the number of real peaks)
    for i in 1:length(x_troughs)

        # Check if trough is less than 50% of the smaller adjacent peaks
        if y_troughs[i]*2 < min(y_peaks[i], y_peaks[i+1])

            push!(x_peaks_significant, x_peaks[i], x_peaks[i+1])
            push!(y_peaks_significant, y_peaks[i], y_peaks[i+1])
            push!(x_troughs_significant, x_troughs[i])
            push!(y_troughs_significant, y_troughs[i]) 

        end
        
    end

    return(x_peaks_significant, y_peaks_significant, x_troughs_significant, y_troughs_significant)

end    # eliminate_minor_peaks()


function filter_spectral_values(peak_values, peak_freqs)
################################################
    # remove nearby lesser peaks
    
    filtered_peak_values = Float64[]  # Initialize an empty array for filtered values
    filtered_peak_freqs = Float64[]   # Initialize an empty array for filtered frequencies
    n = length(peak_values)

    for i in 1:n-1
        freq_diff = peak_freqs[i+1] - peak_freqs[i]
        if freq_diff > 0.02
            push!(filtered_peak_values, peak_values[i])
            push!(filtered_peak_freqs, peak_freqs[i])
        end
    end

    # Add the last value (since there's no next value to compare)
    push!(filtered_peak_values, peak_values[n])
    push!(filtered_peak_freqs, peak_freqs[n])

    return(filtered_peak_values, filtered_peak_freqs)

end    # filter_spectral_values()


function find_peaks_and_valleys(t, y, peak_vals)
################################################
    
"""
    Calls: eliminate_minor_peaks()
    Called by: Main
"""    

    freqs = t
    spectrum = y

    maxima = findmaxima(spectrum)
    pidx = maxima[1]
    peak_values = maxima[2]
    peak_freqs = freqs[pidx]
    
#    peak_freqs,peak_values =  filter_spectral_values(peak_values, peak_freqs)

    # Set a minimum value (20% of fp value) and only accept peaks above this
    max_peak = maximum(peak_vals)
    cutoff = max_peak * 0.02
##    println(cutoff)

    # Only keep peaks above cutoff
    filtered_peak_indices = findall(x -> x >= cutoff, peak_values)
    filtered_peak_values = peak_values[filtered_peak_indices]
    filtered_peak_freqs = peak_freqs[filtered_peak_indices]

    # Merge peaks
    merged_freqs = []
    merged_values = []
    tolerance = 0.05
    for i in 1:length(filtered_peak_freqs)
        if isempty(merged_freqs) || minimum(abs.(filtered_peak_freqs[i] .- merged_freqs)) > tolerance
            push!(merged_freqs, filtered_peak_freqs[i])
            push!(merged_values, filtered_peak_values[i])
        else
            for j in 1:length(merged_freqs)
                if abs(filtered_peak_freqs[i] - merged_freqs[j]) <= tolerance
                    if filtered_peak_values[i] > merged_values[j]
                        merged_freqs[j] = filtered_peak_freqs[i]
                        merged_values[j] = filtered_peak_values[i]
                    end
                end
            end
        end
    end

    # Finding valleys/separation points between merged peaks
    separation_freqs = []
    separation_values = []
    for i in 2:length(merged_freqs)
        # Check interval between two successive peaks
        interval_idx = findall(x -> x>=merged_freqs[i-1] && x<=merged_freqs[i], freqs) 

        if !isempty(interval_idx)
            interval_spectrum = spectrum[interval_idx]
            minima = findminima(interval_spectrum)

            # find minimum value in valleys in the separation point
            if !isempty(minima[1])
                # findminima returns all minima in interval. We only want the lowest.
                min_val_index = argmin(minima[2])
                valley_idx = interval_idx[minima[1][min_val_index]]
                push!(separation_freqs, freqs[valley_idx])
                push!(separation_values, minima[2][min_val_index])
            end
        end
    end
    
    merged_freqs, merged_values, separation_freqs, separation_values = eliminate_minor_peaks(merged_freqs, merged_values, separation_freqs, separation_values)

    return(merged_freqs, merged_values, separation_freqs, separation_values, cutoff)
    
end    # find_peaks_and_valleys()

# use a nearest neighbour approach to locate the closest point
function nearest_neighbor(target, points)
#########################################
    
"""
    Calls: Nil
    Called by: get_wave_types()
"""
    
    distances = [norm(target - point) for point in points]
    nearest_index = argmin(distances)
    
    return points[nearest_index]
    
end


function get_wave_types(start, next, t, y, dir, p1, previous_wave_type, separation_point, sea_swell_transition)
##############################################
    
"""
    Calls:  calc_moments()
            nearest_neighbor()
    Called by: Main
"""

    t_ = t[start:next]
    y_ = y[start:next]
    dir_ = dir[start:next]

    p1 = Plots.plot!(t_,y_, lw=:2, fillrange=0, fillalpha=0.1, label="$(t[start]) to $(t[next])Hz.\n")
    p1 = Plots.vline!([t[next]], label="")
    
    peak = argmax(y_)
    fₚ = t_[peak]
    Tₚ = 1/fₚ    
    m₋₁, m₀, m₁, m₂, m₄ = calc_moments(t_,y_)
    hₘ₀ = 4√m₀
    Hᵣₘₛ = √(8m₀)
    T₀₂ = √(m₀/m₂)    # mean wave period
    T₀₁ = m₀/m₁       # significant wave period
    T₋₁₀ = m₋₁/m₀     # mean energy period Tₑ
    Tc = √(m₂/m₄)
    α = 0.0081
    
    # get peak direction and mean direction (weighted by spectral energy values)
    peak_direction = dir_[peak]
    weighted_mean_direction = mean(dir_, Weights(y_))

    # Calculate representative P-M spectra
    Sf = [α*g^2 * (2π)^-4 * ff^-5 * exp(-1.25 * (ff/fₚ)^-4) for ff ∈ t_]
    p1 = Plots.plot!(t_, Sf, lc=:grey, lw=:2, ls=:dot, label="")

    # determine whether this part of wave record is sea or swell
    # Refer https://www.researchgate.net/publication/249605099_Spectral_Partitioning_and_Identification_of_Wind_Sea_and_Swell
    ratio = y_[peak] / Sf[peak]
    wave_type = (ratio < 1 ? "Ocean Swell" : "Wind Sea")
##    println(ratio)
    # test whether separation between sea and swell has been found
    if (previous_wave_type == "Ocean Swell") && (wave_type == "Wind Sea")
    # Separation between swell and sea has been located
        separation_point = start
        sea_swell_transition = separation_point
##        println(sea_swell_transition)
    end

    previous_wave_type = wave_type
       
    # get direction string
    compass = dir_string[findfirst(x->x==nearest_neighbor(peak_direction, dir_brgs), dir_brgs)]
       
##    @printf("Hₘ₀ = %5.2fm; Hᵣₘₛ = %5.2fm; T₋₁₀ = %5.2fs; T₀₁ = %5.2fs; T₀₂ = %5.2fs; Tc = %5.2fs; Tₚ = %5.2fs; Peak Dirn = %6.2fᵒ; Mean Dirn = %6.2fᵒ Wave type = %s from %s\n",
##            hₘ₀, Hᵣₘₛ, T₋₁₀, T₀₁, T₀₂, Tc, Tₚ, peak_direction, weighted_mean_direction, wave_type, compass)
##    flush(stdout)
    
    return(previous_wave_type, separation_point, sea_swell_transition)
        
    end    # get_wave_types()


# find closest frequency in arr to target value
function search_nearest_index(arr, target)
##########################################

differences = abs.(arr .- target)
nearest_index = argmin(differences)

return(nearest_index)

end    # search_nearest_index()


function calc_model(model, frequency, spectra)
##############################################

    # Initial guess for the parameter
    p0 = [1.0]  # Adjust as needed

    fit_result = curve_fit(model, frequency, spectra, p0)

    # Extract fitted parameters for f^-x curve
    amplitude = fit_result.param[1]

    spectral_fit = model(frequency, fit_result.param)

    return(spectral_fit)

end    # calc_model()


