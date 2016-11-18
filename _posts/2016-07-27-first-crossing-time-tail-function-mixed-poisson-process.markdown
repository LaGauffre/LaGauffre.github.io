---
layout: post
title: When a mixed Poisson process is chasing an upper boundary
---

###Introduction

In a boundary crossing problem, we have

1. a stochastic process that models the random evolution of a quantity over time,
2. an upper boundary, that is a deterministic function of time, usually non-decreasing.

Many problems in applied probabilty comes down to the study of the distribution of the stopping time corresponding to the rendez-vous of the stochastic process and the upper boundary. Typical field of aplications are sequential analysis, reliability, queueing theory, and the one being at the center of the present post - Risk theory.

In risk theory, a non-life insurance company, operating for instance on the car insurance market, is assumed to be able to follow the evolution of its financial reserves continuously in time. As time goes by, some unfortunate policy-holders are facing hazardous events, a car accident say, causing corporal or material damages. The role of the insurance company is to offer a financial compensation in case of the occurence of a claim. Up to a time horizon $$t\geq0$$, the number of claims is modeled throuhg a counting process $$\{N(t)\text{ ; }t\geq0\}$$. The classical assumption consists in opting for a homogeneous Poisson process. The losses, from the insurance company point of view, associated to each claim form a sequence $$\{U_k\text{ ; }k\geq1\}$$ of **i.i.d.**, non-negative random variables. The expenses of the insurance company at time $$t$$ are of amounts $$\left\{X(t)=\sum_{k=1}^{N(t)}U_k \text{ ; }t\geq0\right\}$$. To benefit from the services of the insurance company, the customers pay a premium. The large number of policy-holders in the portfolios enables to approximate the premium flow as a linear function of time with a slope $$c\geq0$$. The financial health and stability of the insurance company coincides with a wealth remaining positive over time. A negative financial reserve is synonymous of bankruptcy. In order to avoid an early ruin situation, the insurance company holds an initial capital of amount $$u\geq0$$. Adding up all the aforementionned elements yields the risk reserve process $$\{R(t)\text{ ; }t\geq0\}$$, defined as

$$
R(t)=u+ct-\sum_{k=1}^{N(t)}U_k.
$$

Herafter a visualization of the process
$$

\begin{tikzpicture}
  %Origin and axis
  \coordinate (O) at (0,0);
  %Former axis
  \draw[dotted,black,->] (0,-0.5) -- (0,9) coordinate[label = {right:$s$}] (xmax);
  \draw[dotted,black,->] (3,0) -- (-5,0) coordinate[label = {below:$y$}] (ymax);
  \draw[black,->] (3,1.5) -- (-5,1.5) coordinate[label = {below:$t$}] (xmax);
  \draw[black,->] (0,0) -- (0,8.5) coordinate[label = {right:$x$}] (ymax);
  %upper linear boundary
  \draw[thick,red] (2,0) -- (-5,9.33);
  \draw (2,0) node[red,below] {$-v$} node{\color{red}$\bullet$};
  \draw (0,2.66) node[red,above right] {$u=\frac{v}{d}-\Delta S_1$} node{\color{red}$\bullet$};
  %Stochastic process trajectory
  \draw[very thick,dotted,blue] (0,0) -- (0,1.5) node[pos=0.5, left] {$\Delta S_1$} ;
  \draw[very thick,blue] (0,1.5) -- (-1.5,1.5) node[pos=0.5, above] {$Y_1$};
  \draw[very thick,dashed,blue] (-1.5,1.5) -- (-1.5,2.5) node[pos=0.5, left] {$\Delta S_2$};
  \draw[very thick,blue] (-1.5,2.5) -- (-3,2.5) node[pos=0.5, above] {$Y_2$};
  \draw[very thick,dashed,blue] (-3,2.5) -- (-3,5) node[pos=0.5, left] {$\Delta S_3$};
  \draw[very thick,blue] (-3,5) -- (-3.75,5) node[pos=0.5, above] {$Y_3$};
  \draw[very thick,dashed,blue] (-3.75,5) -- (-3.75,9) node[pos=0.5, left] {$\Delta S_4$};
  % %Jump Times
  % \draw (0,1.5) node[blue,left] {$S_1$} node{ \color{blue}$\bullet$};
  \draw (0,2.5) node[blue,right] {\hspace{0.3cm}$S_2-\Delta S_1$} node{ \color{blue}$\bullet$};
  \draw (0,5) node[blue,right] {$S_3-\Delta S_1$} node{ \color{blue}$\bullet$};
  Aggregated Capital gains
  \draw (-1.5,1.5) node[blue,below] {$\nu_1$} node{\color{blue}$\vert$};
  \draw (-3,1.5) node[blue,below right] {$\nu_2$} node{ \color{blue}$\vert$};
  % \draw (-3.75,1.5) node[blue,below] {$$} node{ \color{blue}$\vert$};
  %Ruin time in the insurance risk model = First-crossing time
  \draw (-3.75,1.5) node[black,below] {$\tau_u=\nu_3$} node{\color{black}$\times$};
  \draw[dotted,black] (-3.75,1.5) -- (-3.75,7.66);
  \draw[dotted,black] (-3.75,7.66) -- (0,7.66);
  %Ruin time in the dual risk model = First-crossing time
  \draw (0,7.66) node[black,right] {$\sigma_v-\Delta S_{1}=\frac{\tau_u}{d}+u$} node{\color{black}$\times$};
\end{tikzpicture}
$$

We focus our interest on the distribution of the first-crossing time $$\tau_\beta=\{t\geq0\text{ ; }X(t)\geq h_{\beta}(t)\}$$, when the stochastic process passes through the boundary. In our setting, the stochastic process is of jumping type. At random times $$\{T_n\text{ ; }n\geq1\}$$, it makes jumps of random heights $$\{U_n\text{ ; }n\geq1\}$$ in order to get closer to the moving boundary, we assume additionally that the process starts at the origin of the axis, i.e. $$X(0)=0$$. The figure bellow provides a visualization of the crossing problem.
![FirstCrossingTime](/Photos/FirstCrossingTimeBlogPost/FirstCrossingTimeBlogPost.png "The first-crossing time of a stochastic process and an upper moving barrier")
One may note that
the crossing occurs at $$\tau_{\beta}=T_5$$, that is to say while the process is jumping.

The study of the distribution of stopping time such as $$\tau_\beta$$ is of interest in many fields of applied probability. I only describe here the application in risk theory. A non-life insurance company, operating for instance on the car insurance market, is assumed to be able to follow the evolution of its financial reserve in continuous time. Up to some time-horizon $$t\geq0$$, some policyholders might have been facing an unfortunate event, say a car accident. The number of reported car accident until time $$t$$ is a counting process $$\{N(t)\text{ ; }t\geq0\}$$, those happen sequentially at random time $$\{T_{n}\text{ ; }n\geq0\}$$ that corresponds to the aforementionned jump times. A car accident induces a loss for the insurance company that must compensate for the corporal or material damages linked to the accident. The amount of the loss is correlated to the magnitude of the event, and modeled through a sequence of non-negative random variables $$\{U_n\text{ ; }n\geq0\}$$. The liability of the insurer is therefore given by

$$X(t)=\sum_{k=1}^{n}U_{n},$$

at time $$t\geq0$$. The insurance company holds an initial capital $$\beta\geq0$$ and is collecting premiums from its customers, the financial reserve increases as $$h(t)$$. By adding up all this elements, we define the risk reserve process $$\{R(t)\text{ ; }t\geq0\}$$ with

$$
R(t)=h_{\beta}(t)-X(t)
$$

for $$t\geq0$$, which gives the level of the financial reserve of the insurance company. Risk theory focuses on the computation of the ruin probability defined as

$$
\psi(\beta,t)=\mathbb{P}(\tau_{\beta}\leq t)
$$

where $$\tau_{\beta}=\inf\{t\geq0\text{ ; }R(t)\leq 0\}$$ is the ruin time, i.e. the first instant at which the financial reserve becomes negative. The probability $$\psi(\beta,t)$$ is the probability that the insurance company goes bankrupt before time $$t$$, given an initial held capital of amount $$\beta$$. The time to ruin can be rewritten as $$\tau_{\beta}=\inf\{t\geq0\text{ ; }X(t)\geq h_{\beta}(t)\}$$. Thus the ruin problem reduces to a boundary crossing problem involving the stochastic process $$\{t\geq0\text{ ; }X(t)\}$$ and the upper barrier $$h_{\beta}(t)$$. The ruin probability is computed for risk management purposes. The premium income and he initial capital are calibrated so that the probability of ruin is low.

In this post, a boundary crossing problem is solved in a rather simple situation where

- the claim sizes are all equal to $$1$$,
- the premium are collected linearly in time at a rate $$c\geq0$$,
- no initial capital is held by the insurance company,
- the claim frequency is modeled through a mixed Poisson process.

The ruin probability has been derived a long time ago in this setting for the homogeneous Poisson process. We propose in this post an alternative proof allowing to derive the exact same formula for the mixed Poisson process. The next section gives some background on mixed Poisson processes and also the mathematical tool used to get the result.

### On the mixed Poisson process and Appell polynomials

### A ballot-type result related to the distribution of the first-crossing time
