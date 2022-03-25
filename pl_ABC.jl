### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 5cb8f806-9c04-11ec-0ee7-e580f31402c7
md"# Approximate Bayesian computation

For a model $M(\phi )$ with a particular property: given
that we have observations $D_{obs}$ , it must be possible to calculate
$p(D_{obs} |\phi)$, the probability of obtaining the observed data, for each
possible model parameterization $\phi$. We will use the term likelihood
synonymously with $p(D_{obs} |\phi)$. On the basis of this
probability, one can derive statistical methods for parameter estimation, model selection and uncertainty analysis.
"

# ╔═╡ 6bc01550-09ad-4a1f-8542-df46b2bde366
md"For simple stochastic processes, the probability $p(D_{obs} |\phi)$ can be
calculated directly. One refers to this property by saying that the
process has a *tractable likelihood*. Practically all statistical models that are used in ecology and biology make assumptions that result in tractable
likelihoods. Regression models, for example, typically assume that the
data were observed with independent observation errors that follow a fixed, specified distribution. As the errors are independent,  $p(D_{obs} |\phi)$
simply separates into the probabilities of obtaining the individual data
points, which greatly simplify the calculation. We
will therefore use the term *statistical model* as a synonym for a stochastic
model with a *tractable likelihood*."

# ╔═╡ 83a2c63d-e3ea-4506-b16a-235809b11d48
md"A stochastic simulation is an algorithm that creates samples from a
potentially complex stochastic process by explicitly sampling from all
its sub-processes (Figs 1 and 2). This sampling allows researchers to
model stochastic ecological processes exactly as they are known or conjectured without having to concentrate on the mathematical
tractability of the conditional probabilities that would need to be
calculated to keep the probability $p(D_{obs} |\phi)$ tractable. Stochastic
simulation models are therefore especially useful for describing
processes where many entities develop and interact stochastically,"

# ╔═╡ cc20f816-65ef-4d7e-a555-aa2f8f1053f0
md"Hence, the crucial difference between a typical statistical model and
a stochastic simulation model is not the model structure as such. Both
are representations of a stochastic process. However, while typical
statistical models allow the calculation of $p(D_{obs} |\phi)$ directly (tractable
likelihood), stochastic simulation models produce random draws $D_{sim}$ 
from the stochastic process by means of simulation."

# ╔═╡ db34afd2-c8b1-4d45-8ba7-9fb80324e80b
md"Despite different origins and little apparent overlap, most of these
methods use the same three essential steps:

>(1) The dimensionality of the data is reduced by calculating summary statistics of observed and simulated data.

>(2) Based on these summary statistics,  $p(D_{obs} |\phi)$, the likelihood of obtaining the observed data $D_{obs}$ from the model $M$ with parameters $\phi$, is approximated.

> (3) For the computationally intensive task of estimating the shape of the approximated likelihood as a function of the model parameters, state-of-the-art sampling and optimization techniques are applied."

# ╔═╡ 6c8fcb59-f66a-4b3e-9abf-492771ab7643
md"These steps allow the linkage of stochastic simulation models to
well-established statistical theory and therefore provide a general
framework for parameter estimation, model selection and uncertainty
estimation by comparison of model output and data (inverse
modelling)."

# ╔═╡ a8f7308e-2234-4455-9ffe-b3e39b8ccee7
md"The first step for comparing stochastic simulation models with
observations is to reduce the dimensionality of simulated and
observed data.

Within this article, we use the term *summary statistic* for such an
aggregation of model output and observed data."

# ╔═╡ 797e5ba1-15e6-473d-b0cf-42d288582d9d
md"## Sufficiency of summary statistics

The idea that data may often be reduced without losing information
for the purpose of statistical inference is known as sufficiency: a
summary statistic (or a set of summary statistics) is sufficient if it
produces an aggregation of the data that contains the same
information as the original data for the purpose of parameter
estimation or model selection of a model or a set of models.
"


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╠═5cb8f806-9c04-11ec-0ee7-e580f31402c7
# ╟─6bc01550-09ad-4a1f-8542-df46b2bde366
# ╟─83a2c63d-e3ea-4506-b16a-235809b11d48
# ╟─cc20f816-65ef-4d7e-a555-aa2f8f1053f0
# ╟─db34afd2-c8b1-4d45-8ba7-9fb80324e80b
# ╟─6c8fcb59-f66a-4b3e-9abf-492771ab7643
# ╟─a8f7308e-2234-4455-9ffe-b3e39b8ccee7
# ╠═797e5ba1-15e6-473d-b0cf-42d288582d9d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
