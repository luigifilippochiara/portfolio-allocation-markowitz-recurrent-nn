# Portfolio asset allocation strategies: from Markowitz to RNNs
[![License][licence-badge]](/LICENSE)

Research project to explore different approaches for optimal portfolio allocation starting from 18 EU bond indices and a benchmark. Project developed using Matlab and Python.

![Portfolio allocation](https://github.com/luigifilippochiara/portfolio-allocation-markowitz-recurrent-nn/blob/main/preview.jpg?raw=true "Portfolio allocation")

## Input data & techniques used in the project
Raw data is all returns, all maturities bond indices prices of 18 EU countries from 1998 to 2018, downloaded from Eikon. Exchange rates were downloaded as well. Techniques used range from: weight budgeting in Markowitz framework, risk budgeting, constant correlation models and recurrent neural networks.

## Project structure
* [main.m](https://github.com/luigifilippochiara/portfolio-allocation-markowitz-recurrent-nn/blob/main/main.m) contains the code to extract the data from the original .xlsx file
* [a001_DataAnalysis.m](https://github.com/luigifilippochiara/portfolio-allocation-markowitz-recurrent-nn/blob/main/a001_DataAnalysis.m) performs the preliminary data exploration and analysis
* [a002_a_StrategicAssetAllocation.m](https://github.com/luigifilippochiara/portfolio-allocation-markowitz-recurrent-nn/blob/main/a002_a_StrategicAssetAllocation.m) performs some initial analysis on the unconstrained and constrained efficient frontier over two different periods of 5 years
* [a002_b_Forecast.m](https://github.com/luigifilippochiara/portfolio-allocation-markowitz-recurrent-nn/blob/main/a002_b_Forecast.m) replicates the possible evolution of the equally weighted and benchmark portfolio over a period of 5 years
* `a003_...m` and `a004_...m` files explores different advanced techniques for portfolio allocation, computing the evolution of the assets weights over time, cumulative returns and overall ranking of the various strategies
* The folder `/RNN` contains the python files related to the recurrent neural network used in `a004_a_Advanced.m`

## Additional info on using the project files
Most of the Matlab `a00x_...m` files should be self contained and use `.mat` files to gather the data needed to perform the analysis contained in the respective files.

## License
Use as you wish. This project is licensed under the MIT License.


[licence-badge]: https://img.shields.io/npm/l/express.svg
