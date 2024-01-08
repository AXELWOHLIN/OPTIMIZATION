import numpy as np
from scipy.optimize import minimize
from cirpdf import cirpdf
import pandas as pd

#get values from pdf -> log values -> return negative log values
#so that minimization algorithm can be used. Check for 0 values
#since log function does not accept these.
def objective_function(params, data, times):
    a, b, sigma = params
    pdf_vals = cirpdf(data[1:], times[1:], a, b, sigma)
    if np.any(pdf_vals <= 0):
        return 1000 
    log_likelihood = np.sum(np.log(pdf_vals))
    return -log_likelihood 

#using pandas to read dataset
df = pd.read_csv('CIRDataSet.csv')
times = df['t'].values
data = df['data'].values

results = []
for a in [0.001,0.01,0.1,0.3,0.5,1,3,5,10,15,20,99]:
    a = a
    b = np.mean(data)
    sigma = np.std(data)
    print(a,b,sigma)
    #give initial start values for params
    #initial_guess = [0.001, 0.2, 0.02]
    initial_guess = [a,b*2,sigma*2]


    bounds = [(0.00001,100), (0.000001, 100), (0.000001, 100)]

    #use scipy minimize 
    result = minimize(objective_function, initial_guess, args=(data, times), bounds=bounds)
    optimized_params = result.x
    final_vals = cirpdf(data, times, optimized_params[0], optimized_params[1], optimized_params[2])
    negative_L_value = -np.sum(np.log(final_vals))

    results.append({
        "a": optimized_params[0],
        "b": optimized_params[1],
        "sigma": optimized_params[2],
        "Negative L-value": negative_L_value
    })

#results to pandas datafram
results_df = pd.DataFrame(results)

#show results
print(results_df)

#save results
results_df.to_csv("optimization_results.csv", index=False)
"""
    #print out results
    print(result.success)
    print(result.message)
    print(f"Optimal parameters: {optimized_params}")
    final_vals = cirpdf(data,times,optimized_params[0],optimized_params[1],optimized_params[2])
    print(f"Negative L-value: {-np.sum(np.log(final_vals)) }")

"""
        
      