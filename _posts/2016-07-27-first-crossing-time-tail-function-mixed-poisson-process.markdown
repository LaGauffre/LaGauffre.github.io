---
layout: post
title: When a mixed Poisson process is chasing an upper boundary (Work in progress!!!)
---

###Introduction

In a boundary crossing problem, we have

1. a stochastic process $$\{X(t)\text{ ; }t\geq0\}$$, that models the random evolution of some quantity over time,
2. a boundary, that is a deterministic function of time.

In this post, we assume that the stochastic process starts at the origin of the axis, so $$X(0)=0$$. The boundary is an upper boundary defined as $$h_{\beta}(t)=h(t)+\beta$$, where $$\beta\geq0$$ and $$h(t)$$ is a non-decreasing function of time. We are interested in the distribution of the first-crossing time $$\tau_\beta=\{t\geq0\text{ ; }X(t)\geq h_{\beta}(t)\}$$, at which the stochastic process passes through the boundary. In our setting, the stochastic process is of the jumping type. At random times $$\{T_n\text{ ; }n\geq1\}$$, $$X(t)$$ jumps. The heights of the jumps form a sequence $$\{U_n\text{ ; }n\geq1\}$$ of **i.i.d.** non-negative random variables. The figure bellow provides a visualization of the crossing problem.
![FirstCrossingTime](/Photos/FirstCrossingTimeBlogPost/FirstCrossingTimeBlogPost.png "The first-crossing time of a stochastic process and an upper moving barrier")
The crossing occurs at $$\tau_{\beta}=T_5$$, while the process is jumping.

The study of the distribution of stopping times such as $$\tau_\beta$$ is of interest in many fields of applied probability. I only describe here one application within ruin theory (Insurance company risk management). A non-life insurance company, operating for instance on the car insurance market, is assumed to be able to follow the evolution of its financial reserves in continuous time. Up to some time-horizon $$t\geq0$$, some policy-holders might have been facing an unfortunate event, say a car accident. The number of reported car accident until time $$t$$ is $$N(t)$$ ($$\{N(t)\text{ ; }t\geq0\}$$ is a counting process). The $$n^{\text{th}}$$ car accident happens at time $$T_{n}$$ and corresponds to the aforementionned jump time. A car accident induces a loss for the insurance company that must compensate for the corporal or material damages linked to the accident. The amount of the loss is correlated to the magnitude of the event, the loss associated to the $$n^{\text{th}}$$ accident is a non-negative random variable denoted by $$U_{n}$$. The claim sizes form a sequence $$\{U_{n}\text{ ; }n\geq1\}$$ of **i.i.d.** non-negative random variables that one can identify to the above jump heights. The liability of the insurer is therefore given by

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

### On the mixed Poisson process and Appell polynomials

### A ballot-type result related to the distribution of the first-crossing time
