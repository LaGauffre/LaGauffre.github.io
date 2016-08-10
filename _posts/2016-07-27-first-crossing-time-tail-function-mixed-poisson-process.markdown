---
layout: post
title: When a mixed Poisson process is chasing an upper boundary
---

###Introduction

In a boundary crossing problem, there are two main characters,

1. a stochastic process $$\{X(t)\text{ ; }t\geq0\}$$, that models the random evolution of a quantity over time,
2. an upper boundary, that is a deterministic function of time $$h_\beta(t)=h(t)+\beta$$, where $$h(t)$$ is non-decreasing and $$\beta\geq0$$.

We focus our interest on the distribution of the first-crossing time $$\tau_\beta=\{t\geq0\text{ ; }X(t)\geq h_{\beta}(t)\}$$, when the stochastic process passes through the boundary. In our setting, the stochastic process is of the jumping type. At random times $$\{T_n\text{ ; }n\geq1\}$$, it makes jumps of random heights $$\{U_n\text{ ; }n\geq1\}$$ in order to get closer to the moving boundary, we assume additionally that the process starts at the origin of the axis, i.e. $$X(0)=0$$. The figure bellow provides a visualization of the crossing problem.
![FirstCrossingTime](/Photos/FirstCrossingTimeBlogPost/FirstCrossingTimeBlogPost.png "The first-crossing time of a stochastic process and an upper moving barrier")
One may note that
the crossing occurs at $$\tau_{\beta}=T_5$$, that is to say while the process is jumping.

The study of the distribution of stopping time such as $$\tau_\beta$$ is of interest in many fields of applied probability. I only describe here the application in ruin theory (Insurance company risk management). A non-life insurance company, operating for instance on the car insurance market, is assumed to be able to follow the evolution of its financial reserve in continuous time. Up to some time-horizon $$t\geq0$$, some policy-holders might have been facing an unfortunate event, say a car accident. The number of reported car accident until time $$t$$ is a counting process $$\{N(t)\text{ ; }t\geq0\}$$, those happen sequentially at random time $$\{T_{n}\text{ ; }n\geq0\}$$ that corresponds to the aforementionned jump times. A car accident induces a loss for the insurance company that must compensate for the corporal or material damages linked to the accident. This amount of the loss is correlated to the magnitude of the event, and modeled through a sequence of non-negative random variables $$\{U_n\text{ ; }n\geq0\}$$. The liability of the insurer is therefore given by

$$X(t)=\sum_{k=1}^{n}U_{n},$$

at time $$t\geq0$$. The insurance company holds an initial capital $$\beta\geq0$$ and is collecting premiums from its customers, the financial reserve as the non-decreasing function $$h(t)$$. Adding up all this elements allows to define the risk reserve process $$\{R(t)\text{ ; }t\geq0\}$$ as

$$
R(t)=h_{\beta}(t)-X(t)
$$

for $$t\geq0$$, which gives the level of the financial reserve of the insurance company. Ruin theory focuses on the computation of the ruin probability defined as

$$
\psi(\beta,t)=\mathbb{P}(\tau_{\beta}\leq t)
$$

where $$\tau_{\beta}=\inf\{t\geq0\text{ ; }R(t)\leq 0\}$$ is the ruin time, i.e. the first instant at which the financial reserve becomes negative. The probability $$\psi(\beta,t)$$ is the probability that the insurance company goes bankrupt before time $$t$$, given an initial held capital of amount $$\beta$$. The time to ruin can be rewritten as $$\tau_{\beta}=\inf\{t\geq0\text{ ; }X(t)\geq h_{\beta}(t)\}$$. Thus the ruin problem reduces to a boundary crossing problem involving the stochastic process $$\{t\geq0\text{ ; }X(t)\}$$ and the upper barrier $$h_{\beta}(t)$$. The ruin probability is computed for risk management purposes. The premium income and he initial capital are calibrated so that the probability of ruin is low.

In this post, a boundary crossing problem is given in a rather simple situation where

- the claim sizes are all equal to $$1$$,
- the premium are collected linearly in time at a rate $$c\geq0$$,
- no initial capital is held by the insurance company,
- the claim frequency is modeled through a mixed Poisson process.

The ruin probability has been derived a long time ago for in this very setting for the homogeneous Poisson process. We propose in this post an alternative proof allowing to derive the exact same formula for the mixed Poisson process. The next section gives some background on mixed Poisson processes and also the mathematical tool used to get the result.

###On the mixed Poisson process and Appell polynomials

###A ballot-type result related to the distribution of the first-crossing time
