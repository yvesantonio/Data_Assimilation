# Data-Assimilation

## Introduction

Data assimilation is a mathematical discipline that seeks to optimally combine theory (usually in the form of a numerical model) with observations. There may be a number of different goals sought, for example—to determine the optimal state estimate of a system, to determine initial conditions for a numerical forecast model, to interpolate sparse observation data using (e.g. physical) knowledge of the system being observed, to train numerical model parameters based on observed data. Depending on the goal, different solution methods may be used. Data assimilation is distinguished from other forms of machine learning, image analysis, and statistical methods in that it utilizes a dynamical model of the system being analyzed.

## Data Assimilation Process

Classically, data assimilation has been applied to chaotic dynamical systems that are too difficult to predict using simple extrapolation methods. The cause of this difficulty is that small changes in initial conditions can lead to large changes in prediction accuracy. This is sometimes known as the butterfly effect – the sensitive dependence on initial conditions in which a small change in one state of a deterministic nonlinear system can result in large differences in a later state.

At any update time, data assimilation usually takes a forecast (also known as the first guess, or background information) and applies a correction to the forecast based on a set of observed data and estimated errors that are present in both the observations and the forecast itself. The difference between the forecast and the observations at that time is called the departure or the innovation (as it provides new information to the data assimilation process). A weighting factor is applied to the innovation to determine how much of a correction should be made to the forecast based on the new information from the observations. The best estimate of the state of the system based on the correction to the forecast determined by a weighting factor times the innovation is called the analysis. In one dimension, computing the analysis could be as simple as forming a weighted average of a forecasted and observed value. In multiple dimensions the problem becomes more difficult. Much of the work in data assimilation is focused on adequately estimating the appropriate weighting factor based on intricate knowledge of the errors in the system.

The measurements are usually made of a real-world system, rather than of the model's incomplete representation of that system, and so a special function called the observation operator (usually depicted by h() for a nonlinear operator or H for its linearization) is needed to map the modeled variable to a form that can be directly compared with the observation.

//## Advanced Data Assimilation Algorithms and Methods
