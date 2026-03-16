import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor

# 1. Define the metadata
parameterNames = ["gamma0", "c0", "beta0", "btoc0", "lnw0", 
                  "gamma1", "c1", "beta1", "btoc1", "lnw1", 
                  "gamma2", "c2", "beta2", "btoc2", "lnw2"]

parameterDefaults = [-26, 20, -26, 0.5, 0, 
                     -19, 20, -19, 0.5, 0, 
                     -15, 20, -15, 0.5, 0]

parameterEnables = [True, True, True, False, False, 
                    True, True, True, False, False, 
                    True, True, True, False, False]

parameterMins = [-90, 3.5, -90, 0.1, -3, 
                 -90, 5.5, -90, 0.1, -3, 
                 -90, 7.5, -90, 0.1, -3]

parameterMaxs = [0, 21, 0, 0.9, 3, 
                 0, 33, 0, 0.9, 3, 
                 0, 45, 0, 0.9, 3]

# 2. Load and merge the data
df_params = pd.read_csv('optimizerParameters.txt', sep=';', skipinitialspace=True).dropna(axis=1, how='all')
df_results = pd.read_csv('optimizerResults.txt', sep=';', skipinitialspace=True).dropna(axis=1, how='all')

df_params.columns = df_params.columns.str.strip()
df_results.columns = df_results.columns.str.strip()
df = pd.merge(df_params, df_results, on='caseID')

# 3. Filter for active parameters
active_indices = [i for i, enabled in enumerate(parameterEnables) if enabled]
active_names = [parameterNames[i] for i in active_indices]
active_mins = np.array([parameterMins[i] for i in active_indices])
active_maxs = np.array([parameterMaxs[i] for i in active_indices])

X_train = df[active_names].values
y_train = df['F'].values

# 4. Find the best known point
best_known_idx = np.argmin(y_train)
best_known_params = X_train[best_known_idx]
best_known_F = y_train[best_known_idx]

# 5. Define the Trust Region (Conservative Bounds)
# We restrict the search to +/- 10% of the total parameter range around the best point
trust_region_percent = 0.10 
param_ranges = active_maxs - active_mins

# Calculate new local bounds, ensuring we don't exceed the absolute global bounds
local_mins = np.maximum(active_mins, best_known_params - (param_ranges * trust_region_percent))
local_maxs = np.minimum(active_maxs, best_known_params + (param_ranges * trust_region_percent))

# 6. Train the Surrogate Model
model = RandomForestRegressor(n_estimators=200, random_state=42)
model.fit(X_train, y_train)

# 7. Optimize only within the Trust Region
n_candidates = 50000
candidates = np.random.uniform(low=local_mins, high=local_maxs, size=(n_candidates, len(active_names)))

predicted_F = model.predict(candidates)
best_candidate_idx = np.argmin(predicted_F)
best_active_params = candidates[best_candidate_idx]
best_predicted_F = predicted_F[best_candidate_idx]

# 8. Reconstruct the full parameter set
suggested_params = []
active_counter = 0

for i, enabled in enumerate(parameterEnables):
    if enabled:
        suggested_params.append(best_active_params[active_counter])
        active_counter += 1
    else:
        suggested_params.append(parameterDefaults[i])

# 9. Output
print("=== Conservative Optimization Results ===")
print(f"Anchored to Best Known F: {best_known_F:.4f}")
print(f"Predicted F of New Set:   {best_predicted_F:.4f}\n")
print(f"Search restricted to {trust_region_percent*100}% of parameter ranges.\n")
print("Suggested Next Parameter Set:")

for name, val in zip(parameterNames, suggested_params):
    print(f"{name.ljust(10)} : {val:.4f}")
