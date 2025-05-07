# import Pkg; Pkg.add("Distributions"); Pkg.add("Optim"); Pkg.add("PyPlot"); Pkg.add("CSV"); Pkg.add("DataFrames")
using Distributions, Optim, CSV, DataFrames

#Log-likelihood
function dist_vec(y_vec,x_vec; m=m, c=c, σ=σ, w=w)
    normal_probs = pdf.(Normal(0.0,σ),y_vec .- (m .* x_vec .+ c))
    return sum(log.((w .* normal_probs) .+ (1-w)/((maximum(y_vec) - minimum(y_vec)))))
end

#Logistic parameter transform for bounded weight parameter
w_trans(x) = 1/(1+exp(-x))

#Wrapper for using in Optim, using exp transforms for the std. dev., and a logistic transform for the weight.
function to_maximum(params, y_vec, x_vec)
    return -dist_vec(y_vec, x_vec, m = params[1], c = params[2], σ=exp(params[3]), w = w_trans(params[4]))
end

#Finds ML parameters, and returns just slope and intercept
function estimate_baseline_mc(y_vals,x_vals; verbose = false)
    res = optimize(b -> to_maximum(b, y_vals, x_vals), [1.0,0.0,log(4.0),0.0]);
    verbose && println(res)
    opt_vec = Optim.minimizer(res)
    verbose && println("weight = ", w_trans(opt_vec[4]))
    return opt_vec[1], opt_vec[2] #We only need the regression line for our correction
end


df = CSV.read("../data/mixmod_in.csv", DataFrame);
new_df = deepcopy(df);

x_column = "Bare.bead.2" #This selects the baseline
# Save baseline parameters
o_m_vec = Vector();
o_c_vec = Vector();
mi_vec = Vector();
ma_vec = Vector();

for y_column in names(df)[3:end]

    #Ignore the blanks for fitting:
    nonblank_inds = (df.Sample .!= "Empty");
    x_vals = log2.(df[:,x_column])[nonblank_inds]
    y_vals = log2.(df[:,y_column])[nonblank_inds]
    o_m, o_c = estimate_baseline_mc(y_vals,x_vals)

    #Reloading without Blank bead exclusions for adjusting and plotting
    x_vals = log2.(df[:,x_column])
    y_vals = log2.(df[:,y_column])
    adjusted_y = y_vals .- (o_m.*x_vals .+ o_c)
    new_df[!,y_column*"_corrected_log2"] = adjusted_y

    mi,ma = minimum(x_vals),maximum(x_vals)
    #plot([mi,ma],[o_m*mi+o_c,o_m*ma+o_c],color = "green")
    append!(o_m_vec, o_m);
    append!(o_c_vec, o_c);
    append!(mi_vec, mi);
    append!(ma_vec, ma);
end
CSV.write("../data/WithAdjusted.csv",new_df)

line_df = DataFrame(protein = names(df)[3:end],
                    o_m_val = o_m_vec,
                    o_c_val = o_c_vec,
                    mi_val = mi_vec,
                    ma_val = ma_vec)

CSV.write("../data/line_df.csv", line_df)
