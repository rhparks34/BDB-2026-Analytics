All code used for the 2026 BDB Analytics competition can be found here.  
The Appendix in the attachments section of the submission organizes the code better and associates what each file does. The solveOptimalRunner.m solves the 
bounded-acceleration time-optimal route between constrained endpoints. The energyConstrainedPath.m and the associated dynamics and cost files solve the 
effort-constrained distance minimization problem. The matlab_test2.ipynb calls the MATLAB solvers and runs them against all plays in the dataset.
The efficient_routes.ipynb notebook contains all tests done on the WR Path Efficiency metric, albiet in a very unorganized way (sorry). 

Files containing final metrics/ predictions:
- preds_matlab.joblib contains routes of all players calculated by solveOptimalRunner.m
- efficiency_leaderboard.joblib contains leaderboard of WR Path Efficiency metric for Tmin > 1.5
- preds_NLP.joblib contains routes of all players calculated by energyConstrainedPath.m
- rmse_NLP.joblib contains RMSE of each player's route for each play in preds_NLP.joblib
- rmse_1000_new_NLP.joblib contains RMSE of each player's route for each play for updated distance hueristic (deceleration allowed)

