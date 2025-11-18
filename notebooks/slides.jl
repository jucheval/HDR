### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ d7a03d41-86af-45fc-895d-b1753d5f2a35
begin
    using PlutoUI
    using PlutoTeachingTools
    using WGLMakie
    using Bonito
    using PointProcesses
    using MeanFieldGraph
    using DrWatson
    using Distributions
    using GLFW
    import Random: seed!, randperm
    import MarkdownLiteral: @markdown
end

# ╔═╡ 4549be08-558f-4773-b9c4-6fee5a60c97d
ChooseDisplayMode()

# ╔═╡ cdc9af88-ddb8-459e-9041-32d9d02a8e42
monitors = GLFW.GetMonitors();

# ╔═╡ a10cde06-6230-4d31-8092-0ad0b4e33a0f
md"Monitor: $(@bind monitor Select(monitors))"

# ╔═╡ 82e5f7b8-6e9c-4ff0-a49a-22daaf8300cf
monitorresolution = (monitor.videomode.width, monitor.videomode.height);

# ╔═╡ 49e5bea3-ee26-4cd1-b0d1-707b2a5ad337
md"# An *animated* journey in the fields of PDE, probabilities and statistics with point processes"

# ╔═╡ 673768b4-e18e-4e68-bd26-315ef9588e08
md"""
### Julien Chevallier
\

##### HDR defense, December 12 2025
"""

# ╔═╡ 174d611c-1184-426f-8a69-4f9cca0cf989
html"""<div style="height: 200px;"></div>"""

# ╔═╡ 34f58220-87c3-43d0-bca7-a7afbd6cd9c6
html"""Powered by &nbsp <img alt="Pluto.jl" src="https://raw.githubusercontent.com/fonsp/Pluto.jl/dd0ead4caa2d29a3a2cfa1196d31e3114782d363/frontend/img/logo_white_contour.svg" height=37 > and <img alt="Makie.jl" src="https://github.com/MakieOrg/Makie.jl/blob/f9b432dd4409e51962a4d331a74db6f3faab06b4/images/logo_with_text.svg?raw=true" height=40>"""

# ╔═╡ cbcef609-0781-4728-8a6a-937f82b67c1a
md"## My scientific research field"

# ╔═╡ dd29b1c5-6bd9-48c6-aaa8-de2e951c44ce
md"""
TODO : refaire les box pour griser RDE.
FIXME:
- 3 mathematical domains evolving around the notions of PP and neuro.
- Manuscript organized by mathematical fields to emphasize this plurality. Slides organized by papers.
- PP/Neuro : Patricia thaught me most of what I know. Explain the color code. Keywords.
- PDE : At the beginning of my PhD I didn't know much about it. Marie and Maria helped me. Keywords.
- Proba : I met Guilherme during postdoc. Mean-field and COUPLING.
- Statistics : Fun fact. I started my career with stats but only came back to this topic recently. Collab with Grenoble.
- Julia enthusiast for 2 years and I am using Pluto.jl and Makie.jl
- Choice: avoid technicalities, favor high level explanations, plots and animations
"""

# ╔═╡ f04d0ada-3fd5-4441-8d07-67dec5da2504
md"## Table of contents"

# ╔═╡ 9bdb9101-ec7d-4065-8a17-27e5676951de
md"""
### 1. [Neurosciences](#Neurosciences-(1/2))
\
	
### 2. [Point processes](#Point-processes)
\

### 3. [Coupling](#Coupling-(1/2))
\

### 4. [Large networks approximation](#Large-networks-approximation)
\

### 5. [Large networks estimation](#A-linear-model-with-inhibition)
"""

# ╔═╡ 8a68c356-c1b7-492b-bdc7-b7ac2dd96a3a
md"""
FIXME: add links to every sections above
- first 3 parts are introductions to central topics in my career
- Then, applied to 2 problems
"""

# ╔═╡ 0ca0d6a0-2bfd-4f29-aeff-a81b87e3be3b
md"## Neurosciences (1/2)"

# ╔═╡ bbbfcb51-81db-4b1d-8382-44c1658658aa
Columns(@markdown("""
### Neuron scale:
- Excitable cell -> action potentials (**spikes**)
- Firing rate is a function of the membrane voltage
- <span id="gray">Spikes travel through axons</span>
- Refractory period (voltage reseting)
		
### Network scale:
- Synaptic connections between neurons
- Populations of <span id="red">excitatory</span>/<span id="blue">inhibitory</span> neurons
- <span id="gray">Excitatory/inhibitory synapses</span>
- Spatial organization (e.g. visual cortex)
- Neural oscillations
- <span id="gray">Neural plasticity (**learning**)</span>
	"""),
    App() do
        tmin = 0.0
        tmax = 10.0
        pad = 0.5
        seed!(1)

        fig = Figure(size=(600, 500))
        sg = SliderGrid(fig[2, 1],
            (label="Number of neurons (n):", range=1:100, startvalue=1),
            tellwidth=false)
        n_obs = sg.sliders[1].value

        h = rand(MultivariatePoissonProcess(ones(100)), tmin, tmax)

        ax = Axis(fig[1, 1], ylabel="neuron", xlabel=L"t", title="Raster plot", subtitle="(unit rate Poisson process)", limits=((tmin - pad, tmax + pad), (1.0 - pad, 1.0 + pad)))

        on(n_obs) do n
            limits!(ax, tmin - pad, tmax + pad, 1.0 - pad, n + pad)
        end

        hideydecorations!(ax, label=false)

        scatter!(ax, event_times(h), event_marks(h), markersize=lift(x -> 10 - 5 * (x - 1) / (100 - 1), n_obs), color=:black)


        DOM.div(fig; style="text-align: center")
    end
)

# ╔═╡ 5a5d3447-fbcc-42a1-92f1-f6bbb6825d3f
md"""FIXME
- 2 scales
- as the voltage of a neuron increases, it charges and then release this energy through a spike
- neurons communicate via spikes
- explain refractory period
- show raster plot -> spike trains
- neuron interactions via synapses
- interplay between E/I populations
- 2 phenomena emmerging at the level of the network : spatial organization, oscillations
"""

# ╔═╡ 6007b1d3-17b0-47e8-ab83-6587a66d17e8
md"## Neurosciences (2/2)"

# ╔═╡ 1dfc7065-05bc-4187-9677-52ecf9a33c1e
md"## Point processes"

# ╔═╡ 4e8e58b5-281f-42ab-b3d1-0725ebb5a34f
md"""
- ``n``-variate point process: ``n`` random sets ``N^i`` of points -> raster plots.
- Simplest model: Poisson point process.
"""

# ╔═╡ 68807259-fbb3-4c10-b7b4-5980a59c57e2
md"""
- Representation of point processes thanks to unit rate Poisson process via:
  - SDE/Thinning
  - Time change
"""

# ╔═╡ 58e967cc-6937-4941-953b-f79c822559c4
md"""FIXME
- spike trains / raster plots
- Poisson generalized
- neuron $i$ excites or inhibits neuron $j$
- main tools I used are the two representations
"""

# ╔═╡ c448e7c4-34c9-43bf-9171-0f9a1915e8fb
md"## Hawkes processes"

# ╔═╡ dfcf884d-63cd-4e99-a81b-8338760fb5d5
@markdown("""
### Intensity of ``n``-variate Hawkes process

```math
		  \\lambda^{n,i}_t = f_i( U^{n,i}_{t-} ) \\quad \\text{ with } \\quad U^{n,i}_t := \\sum_{j=1}^n \\int_{0}^{t-} h_{j\\to i}(t-s) N^{n,j}(\\mathrm{d}s)
```
- Usually, ``f_i \\nearrow`` so that ``{\\color{red} h_{j \\to i} > 0}`` (*resp.* ``{\\color{blue} h_{j \\to i} &lt 0}``) means <span id="red">excitation</span> (*resp.* <span id="blue">inhibition</span>).
- Convenient cases: ``f_i(u) = \\mu + u`` and ``h_{j \\to i}(t) = w e^{-\\alpha t}``.
""")

# ╔═╡ 7da9719f-2617-4eee-86fa-20dba7728dd2
md"""FIXME
- intensity is function of voltage
- counting measure of the spikes
- f linear and h exponential
- move $\mu$, $\alpha$ and $\beta_{22}$
"""

# ╔═╡ a933cf6b-306f-4940-a4ca-bf493376c515
md"## Coupling (1/2)"

# ╔═╡ 07f4c0df-0d77-4794-bbcd-9a9533c47283
Columns(DOM.div(md"""
- Probability distributions: dual of b.c. functions
  - ``X_n \to X`` iff ``\mathbb{E}\, \varphi(X_n) \to \mathbb{E}\, \varphi(X)`` for all ``\varphi``
""",
        md"""
        - Coupling of ``\nu_{1}`` and ``\nu_{2}``:
          - a couple ``({\color{darkblue}X_1}, {\color{orange}X_2})`` with marginals ``{\color{darkblue}\nu_1}`` and ``{\color{orange}\nu_2}``.
          - Examples: independent coupling, quantile coupling $X_2 = F_{X_2}^{-1}(F_{X_1}(X_1))$.
        """,
        md"""
        - Interests: 
          - coupling metrics: Wasserstein, $\mathrm{TV}$, Lévy-Prokhorov
          - numerical illustrations
        """),
    App() do
        fig = Figure(size=(700, 500))

        begin # buttons
            buttongrid = fig[1, 1] = GridLayout(tellwidth=false)
            buttonlabels = ["Independent", "Translation"]
            buttons = buttongrid[1, 1:2] = [Makie.Button(fig, label=buttonlabels[i]) for i in 1:2]

            buttons[1].buttoncolor = Makie.COLOR_ACCENT_DIMMED[]
        end

        begin # observables
            X1 = Observable(rand(Normal(), 4) .+ 0.1)
            X2 = Observable(rand(Normal(), 4) .- 0.1)
        end

        begin # Plot layout parameters
            marksize = 13
            ypad = 0.1
        end

        begin # Plot
            ax = Axis(fig[2, 1], ylabel="X", xlabel="realization")
            limits!(ax, (nothing, nothing), (-3.0 - ypad, 3.0 + ypad))

            scatter!(ax, 1:4, X1, label=L"X_1 \sim \mathcal{N}(\mu_1,1)", markersize=marksize, color=:darkblue)
            scatter!(ax, 1:4, X2, label=L"X_2 \sim \mathcal{N}(\mu_2,1)", markersize=marksize, color=:orange)

            axislegend(ax, position=:ct)
        end

        # Connect buttons to actions
        on(buttons[1].clicks) do _
            X1[] = rand(Normal(), 4) .+ 0.1
            X2[] = rand(Normal(), 4) .- 0.1

            buttons[1].buttoncolor = Makie.COLOR_ACCENT_DIMMED[]
            buttons[2].buttoncolor = RGBf(0.94, 0.94, 0.94)
        end
        on(buttons[2].clicks) do _
            X1[] = rand(Normal(), 4) .+ 0.1
            X2[] = X1[] .- 0.2

            buttons[2].buttoncolor = Makie.COLOR_ACCENT_DIMMED[]
            buttons[1].buttoncolor = RGBf(0.94, 0.94, 0.94)
        end

        DOM.div(fig; style="text-align: center")
    end)

# ╔═╡ 661246fe-aea7-48c9-a49d-aed212b10a32
md"""
!!! note "Wasserstein metric (order p)"
	$W_p(\nu_{1}, \nu_{2}) := \inf_{\pi} \left(  \int_{E \times E} d(x_1,x_2)^p \, \pi(\mathrm{d}x_1,\mathrm{d}x_2) \right)^{1/p}= \inf_{(X_1,X_2)} \mathbb{E}\left[ d(X_1,X_2)^{p} \right]^{1/p},$
"""

# ╔═╡ 3f5dba59-d773-45ac-b7e7-f555374a6340
md"""FIXME
- usually we deal with distributions in a weak sense (tested against functions)
- coupling gives more concrete sense (called strong)
- in terms of r.v., a coupling is a couple of r.v. such that the marginals are prescribed
- examples
"""

# ╔═╡ 9ff8a2e5-4e31-4971-bbd6-fc17d9f1b737
md"## Coupling (2/2)"

# ╔═╡ da872269-de92-4003-a2f5-4bb4a79ea228
md"# Large networks approximation"

# ╔═╡ ff5208f3-7116-4dbe-945c-b0ab66e74f9d
Columns(@markdown("""				  
<h3 id="gray">Age structure <span id="small-gray">(Pakdaman et al., 2010)</span></h3>
				  
<span id="small-gray">
- <em>Mean-field limit of generalized Hawkes processes</em>.</br>
Stoc. Proc. Appli., 2017</br></br>
- <em>Fluctuations for mean-field interacting age-dependent Hawkes processes</em>.</br>
Elec. Journal Prob., 2017</br></br>
- <em>Stimulus sensitivity of a spiking neural network model</em>.</br>
Journal Stat. Phys., 2018
</span>
"""),
    @markdown("""
  ### Spatial structure <span id="small-gray">(Amari, 1977)</span>


  <span id="small">
  - <em>Mean field limits for nonlinear spatially extended Hawkes processes with exponential memory kernels</em></br>
  with A. Duarte, E. Löcherbach, G. Ost </br>
  Stoc. Proc. Appli., 2019</br></br>
  - <em>Fluctuations for Spatially Extended Hawkes Processes</em>.</br>
  with G. Ost</br>
  Stoc. Proc. Appli., 2020</br></br>
  </span>
  <span id="small-gray">
  - <em>Uniform decomposition of probability measures: quantization, clustering and rate of convergence</em>.</br>
  Journal Applied Prob., 2018
  </span>
  """),
    @markdown("""
  ### Multiclass structure <span id="small-gray">(Ditlevsen & Löcherbach, 2017)</span>

  <span id="small">
  - <em>Diffusion approximation of multi-class Hawkes processes: Theoretical and numerical analysis</em></br>
  with A. Melnykova, I. Tubikanec </br>
  Adv. Applied Prob., 2021
  </span>
  """),
    html"""
  		<style>
  	#small {
  		font-size: 20px;
  	}

  	#small-gray {
  		font-size: 20px;
  		color: gray;
  	}
  </style>
  		""",
    widths=[33, 33, 33, 0])

# ╔═╡ 0cfd87f6-3e2f-423a-894a-7a3223b493c8
html"</br>"

# ╔═╡ 40dc64c1-cdfc-49a2-aa79-63aea3e52b3c
html"</br>"

# ╔═╡ dc8f4bd1-77ea-49c9-8dac-54f8f7064d48
md"## Mean-field approximation"

# ╔═╡ 4394c62c-57bf-486a-af2d-cc9c16715e26
md"""FIXME
- green particle is not really influenced by each other particles individually but merely by their statistical distribution
- in the limit, it is influenced by the limit distribution of the particles including itself
- 3 objectives in the manuscript but I preferred 4 today.
- historical survey
- Objective 4 is original
"""

# ╔═╡ eeb75f07-3dbe-4b48-8498-e532f2418ed6
DOM.div(md"### Main objectives",
    md"""
    1. Find the limit equation and associated Fokker-Planck equation [Weiss, 1907]
        - Empirical distribution is replaced by distribution of the limit process
        - Solve McKean-Vlasov non-linearity via the PDE
    """,
    md"""
    2. Prove that the approximation is consistent (rate $n^{-1/2}$) [Sznitman, 1991]
        - Convergence of the empirical distribution is equivalent to convergence of pairs of particles
        - Couple dependent particles with independent limit particles via Thinning
    """,
    md"""
    3. Prove a CLT for the fluctuations [Méléard, 1996]
        - Approximation rates in terms of $W_p$ or $\mathrm{TV}$
        - Tightness criterion in a suitable distributional space
    """,
    md"""
    4. Derive a first-order approximation (rate $n^{-1}$)
        - Fluctuations CLT $\leftrightarrow$ linearization of a non-linear SPDE
    """)

# ╔═╡ a3a74379-384e-4e04-9a69-dbaad58a7cd2
md"## Spatial Hawkes process / Neural Field Equation (1/3)"

# ╔═╡ ee2d968a-a662-4da4-a17c-19839dffbfc2
md"""FIXME
- neurons have positions in $R^d$
- intensity is a function of the voltage
- voltage satisfy an SDE (corresponds to exp. delay function)
- Right : heatmap of the voltage as t and x
- NFE : same dyanmics for the voltage but there is a continuous spatial distribution of neurons. 
- The averaged counting measures is replaced by the firing rate
"""

# ╔═╡ 05f5212c-44a8-40c7-89b8-ad1abf863fc6
md"### Plot parameters"

# ╔═╡ e0e76e79-a3cf-4fc8-841b-67a242da0542
DOM.div(md"""
-  $n=300$
""",
    md"""
    -  $\rho = \operatorname{Unif}([-\pi, \pi])$
    """,
    md"""
    -  $f(u) = (1 + \exp(-(u - u^*) / \kappa))^{-1}$ with $u^* = 1/2$ and $\kappa = 1/20$
    """,
    md"""
    -  $\alpha = 1$
    """,
    md"""
    -  $w(y,x) = 2\pi\cos(y-x)$
    """,
    md"""
    -  $u^{\rm in}(x) = 1.8\cos(x)$
    """)

# ╔═╡ abe2e182-1f82-46ea-8e09-b845e0671878
md"## Spatial Hawkes process / Neural Field Equation (2/3)"

# ╔═╡ 488ffd0d-cfe6-4cbd-9cd5-7a554cc33868
md"""
$\mathcal{F} := \left\{ u \in C(\mathbb{R}_+ \times \mathbb{R}^d),\, \int_{\mathbb{R}^d} \left\{ \int_{0}^t f(u(s,y)) \mathrm{d} s  + \left( \int_{0}^t  f(u(s,y)) \mathrm{d} s \right)^2 \right\} \rho(\mathrm{d} y) < \infty \right\}$
"""

# ╔═╡ dab4e607-41e7-427f-a9a1-a8b43ca6da8e
md"""
!!! info "Objective 1 - CDLO (2018) & CO (2020)"
	Under $\left( \mathcal{A}^{f,w,u^{\rm in}}_{\rm BL} \right)$, assume that $f$ is invertible. Then, (NFE) is weakly well-posed in $\mathcal{F}$.
	
	Assume furthermore that $\rho$ is finite and $\left( \mathcal{A}^{f,w,u^{\rm in}}_{C^\infty_{\rm b}} \right)$ is satisfied. Then, the solution is $C(\mathbb{R}_+, C^\infty(\mathbb{R}^d))$.
"""

# ╔═╡ c4943fad-6567-4b21-a4b2-af9dd2682ba9
Columns(md"3 distributions:",
    md" $\left\{ \lambda_t^{n,i} \text{ at position } x_i \right\}_i$",
    md" $\left\{ f(u(t,x_i)) \text{ at position } x_i \right\}_i$",
    md" $\left\{ f(u(t,x)) \text{ distributed as } \rho \right\}$",
    widths=[16, 28, 28, 28])

# ╔═╡ 6d4ed1ad-37e5-4850-b9a3-8e783b867a16
md"""
!!! tip "Objective 2 - CDLO (2018)"
	Under $\left( \mathcal{A}^{f,w,u^{\rm in}}_{\rm BL} \right)$, assume that $\int e^{\beta \|x\|}\rho(\mathrm{d} x) < \infty$ for some $\beta >0$. Then,

	$\mathbb{E}\left[ \frac{1}{n} \sum_{i=1}^{n} \int_{0}^{t} \left| U^{n,i}_s - u(s,x_i) \right| \mathrm{d} s   \right] \lesssim_{t,f,\alpha,\beta,w,u^{\rm in}}  n^{- 1/2 } +   W_2 ( \rho^n , \rho ).$
"""

# ╔═╡ cda8bffd-793d-4467-bf02-214456982b91
md"## Spatial Hawkes process / Neural Field Equation (3/3)"

# ╔═╡ bf5f5043-da32-4136-b536-0ba4ace19703
md"""
!!! warning "Objective 3 - CO (2020)"
	Under $\left( \mathcal{A}^{f,w,u^{\rm in}}_{C^\infty_{\rm b}} \right)$, assume that $\rho$ is the Lebesgue measure on $[0,1]$. Let $\tilde{U}^{n,i}_t := \sqrt{n} (U^{n,i}_t - u(t,x_i))$ and $\tilde{\nu}^n_t(\mathrm{d}x) := n^{-1} \sum \tilde{U}^{n,i}_t \delta_{x_i}(\mathrm{d}x)$. Then $\tilde{\nu}^n \to {\color{green}\tilde{\nu}}$ which is the solution of $(\mathrm{E})$ in $C(\mathbb{R}_+, C^\infty([0,1])^\prime)$.
"""

# ╔═╡ 0f9ab402-3e7d-4367-a1e3-ead7968fdb03
md"""
$\forall \varphi\in C^\infty([0,1]), \quad \left< \tilde{\nu}_{t},\varphi \right> = e^{-\alpha t} \left< M_{t}, \varphi \right> + \int_{0}^{t} e^{-\alpha(t-s)} \int w(y,x) \varphi(x) f^\prime(u(s,y)) \rho(\mathrm{d}x) \tilde{\nu}_s(\mathrm{d}y)  \mathrm{d} s,
\tag{E}$
where $M$ is a centred Gaussian process in $C^\infty([0,1])^\prime$ whose increments variance is related to $w(y,x)^2 f(u(t,y))$.
"""

# ╔═╡ 34950b24-175e-4f66-a220-fcb7df592f36
md"""
!!! danger "Objective 4 - CO (2020)"
	Furthermore assume that $f$ is lower bounded and let $\hat{u}^n(t,x) := u(t,x) + n^{-1/2} {\color{green}\tilde{\nu}}(t,x)$. 
	We can couple $\hat{u}^n$ and $u^n$ in such a way that for all $t\geq 0$,
	
	$\begin{equation*}
	    \sup_{s\in [0,t], x\in [0,1]} \mathbb{E}\left[ \left| \hat{u}^n(t,x) - u^n(t,x) \right|^2 \right]^{1/2} \lesssim_{t, w, f, \alpha} n^{-1}.
	\end{equation*}$
"""

# ╔═╡ 5012cdc0-2399-4ee0-8819-14c75015d9a2
md"## Multi-class Hawkes process - Ditlevsen & Löcherbach (2017)"

# ╔═╡ bc594f49-72f5-43a9-a5c3-a8694c0b4b89
@markdown("""
!!! warning "Diffusion approximation"
	```math
	\\frac{1}{n_k} \\sum_{i=1}^{n_k} N^k(\\mathrm{d}t) \\quad \\longrightarrow \\quad f_k\\left( U^k_{t-} \\right) \\mathrm{d}t + \\sqrt{\\frac{f_k\\left( U^k_{t-} \\right)}{n_k}} \\mathrm{d} B^k_t
	```
""")

# ╔═╡ 7d3e5eac-cea0-47f7-840e-485dccc715ec
md"""FIXME
- 6 discrete locations
- all to all connection from neurons in pop. 1 to neurons in pop. 2 and so on...
- neurons of the same pop. share the same voltage variable
- Erlang function generalizes the exponential one.
- for each pop. the averaged counting measure is replaced by the firing rate + some diffusion noise
- weak diffusion approximation has been proved by Dit and Loc
"""

# ╔═╡ 2435405c-8428-4af0-8d4e-4ae099cf3f13
md"## Strong diffusion approximation of Poisson process"

# ╔═╡ fc24f7ef-d9c5-4c52-ae27-266850ce49a2
md"""
!!! danger "KMT coupling"
	We can couple a unit rate 1d-Poisson process $\Pi$ and a Brownian motion $B$ in a such a way that
	
	$\sup_{t\geq 0} \frac{|{\color{darkblue}\Pi_t - t} - {\color{orange}B_t}|}{\ln(2 \vee t)} < \infty$
	admits exponential moments.
"""

# ╔═╡ e901002e-e3c7-4d38-be78-c59032d604b6
md"""FIXME
- historical case
-  $B_t$ on a coarse resolution
- I have an animation to explain how it is done if you want
"""

# ╔═╡ b7542278-c465-4096-921d-aee94dc3adbf
md"## Strong diffusion approximation of MC-HP"

# ╔═╡ 0e8bedd2-36e0-4acf-817e-b42a2fd97bbd
@markdown("""
!!! danger "Theorem - CMT (2021) & C (2023)"
	Under ``\\left( \\mathcal{A}^{f_k}_{\\rm BL} \\right)``, we can couple solutions of the <span id="darkblue">Multi-class Hawkes process</span> and the <span id="orange">diffusion approximation</span>:

  	```math
    \\sup_{t\\geq 0} \\frac{\\| {\\color{darkblue}U_t} - {\\color{orange}\\overline{U}_t}\\|_{\\infty}}{e^{\\kappa t}} \\leq K \\frac{\\ln n}{n},
	```

	where ``K`` admits exponential moments.
""")

# ╔═╡ a4551abd-2aa2-47c3-aea4-d23e81af58d8
html"""<div style="height: 100px;"></div>"""

# ╔═╡ bc1554fd-6c3c-4217-9df4-1026d5fe9925
md"### Plot parameters"

# ╔═╡ b38816f9-6cc0-4c16-9d56-cf8f6c2f1ace
DOM.div(md"""
- 2 populations with $n_1 = n_2 = 10$
""",
    md"""
    -  $f_1(x) = 10e^x \mathbf{1}_{(-\infty, \ln 20)}(x) + 400 / (1 + 400e^{-2x}) \mathbf{1}_{[\ln 20, +\infty)}(x)$
    """,
    md"""
    -  $f_2(x) = e^x \mathbf{1}_{(-\infty, \ln 20)}(x) + 40 / (1 + 400e^{-2x}) \mathbf{1}_{[\ln 20, +\infty)}(x)$
    """,
    md"""
    -  $\alpha_1 = \alpha_2 = 1$
    """,
    md"""
    -  $c_1 = -1$, $c_2=1$
    """,
    md"""
    -  $m_1=3$, $m_2=2$
    """,
    md"""
    -  $U^{1,m^\prime}_0 = -2$, $U^{2,m^\prime}_0 = 2$
    """)

# ╔═╡ 89b224b0-fe7d-44d5-a641-bcc4a7a472c3
@markdown"""# Large networks estimation <span id="small-gray">inspired by Delattre & Fournier (2016)</span>"""

# ╔═╡ 533a2768-2c83-483f-9d0d-4dfba4acf8b5
Columns(@markdown("""				  
### Graph density estimation

<span id="small2">
- <em>Inferring the dependence graph density of binary graphical models in high dimension</em>.</br>
with E. Löcherbach and G. Ost</br>
Accepted, AOS
</span>
"""),
    @markdown("""
  ### Community detection

  <span id="small2">
  - <em>Community detection for binary graphical models in high dimension</em></br>
  with G. Ost</br>
  Major revision, AOS
  </span>
  """),
    html"""
  		<style>
  	#small2 {
  		font-size: 25px;
  	}
  </style>
  		""",
    widths=[50, 50, 0])

# ╔═╡ 5abc34bd-5950-4804-8125-be48d3de4449
html"</br>"

# ╔═╡ 1b48f36d-fbd2-497d-a175-43aba5f5a6ff
html"</br>"

# ╔═╡ 51124312-b5e5-4324-92fa-cc00853ea2d6
md"## A linear model with inhibition"

# ╔═╡ 469649cd-a6e5-4d69-8ada-20b566bd3cf1
md"""
- **Random graph:** $\Psi^n = (\Psi^n_{ij})_{1\leq i,j\leq n}$ is a random matrix with i.i.d. entries, $\Psi^n_{ij} \sim \operatorname{Ber}(p)$.
- **Spikes:** $X^{n,i}_t = 1$ if neuron $i$ spikes, $X^{n,i}_t = 0$ otherwise. Conditionally on $\Psi^n = \psi$, the process $X^n$ with values in $\{0,1\}^n$ is a Markov chain such that
$\mathbb{P}_{\psi}\left( X^{n,i}_{t} = 1 | X^n_{t-1} = x \right) = \mu + \gamma \left( \frac{1}{n} \sum_{\color{red}{j\in \mathcal{P}_+}} \psi_{ij} x_j + \frac{1}{n} \sum_{\color{blue}{j\in \mathcal{P}_-}} \psi_{ij} (1-x_j) \right),$
 $\hspace{1em}$ and the coordinates are conditionally independent. The parameters are $\theta := (\mu, \gamma, p)$.
- **Assumption:** $\operatorname{Card}\{\mathcal{P}_+\}/n \to r_+$
"""

# ╔═╡ c60f7a22-703e-4cbe-95d8-09039d080d04
md"""FIXME
- latent random Erdos-Rényi graph
- discrete process in time and values (0,1).
- Once the graph is fixed, the process is a MC
- proportion of excitatory and inhibitory neurons converges to something
- move $\mu$ below 0.5 and leave it low
- move $\gamma$ and leave it low
- move $r_+$
"""

# ╔═╡ c562f965-72e2-4a7b-a649-53ac2d83908a
md"## Graph density estimation"

# ╔═╡ 31f9380c-2888-4d75-b72b-a1bbbf7e3b06
DOM.div(md"""
Let $N^{n,i}_{t}:=\sum_{s=1}^t X^{n,i}_{s}$, $\overline{N}^n_{t}:=n^{-1}\sum_{i=1}^{n}N^{n,i}_{t}$, $\hat{m}^{n,i}_t := N^{n,i}_t/t$ and $\hat{m}^n_t:= \overline{N}^n_{t}/t$. Consider 3 moment estimators:

$\hat{m}^n_t, \ \ \
  \hat{v}^n_t:= \sum_{i=1}^{n}\left(\hat{m}^{n,i}_t - \hat{m}^n_t\right)^2 - \left[ \hat{m}^{n,i}_t - (\hat{m}^{n,i}_t)^2 \right]$ 

$\hat{w}^n_t:=2W_{2\Delta}-W_{\Delta}, 
  \ \text{with} \  
  W_{\Delta}=\frac{n}{t}\sum_{k=1}^{\lfloor t/\Delta \rfloor}\left(\overline{N}^n_{k\Delta}-\overline{N}^n_{(k-1)\Delta}-\Delta\hat{m}^n_t\right)^2.$
""",
    html"""<hr>""",
    md"""
    $\begin{pmatrix}
    \hat{m}^n_t, \hat{v}^n_t,  \hat{w}^n_t
    \end{pmatrix} \xrightarrow[n,t \to \infty]{}
    \begin{pmatrix}
    m,  v,  w
    \end{pmatrix} 
    := \Phi(\theta).$
    The function $\Phi$ is explicit (depends on $r_+$).
    """)

# ╔═╡ 2fc9ed70-0ee5-42ed-b45a-a6acd3e67e91
md"""
!!! danger "Theorem - CLO (2025)"
	If $\Phi(\theta)$ is invertible, define $\hat{\theta}^n_t = \Phi^{-1}(\hat{m}^n_t, \hat{v}^n_t, \hat{w}^n_t)$. Then, for all $\varepsilon\in (0,1)$,

	$\mathbb{P}_\theta\left( \| \hat{\theta}^n_t - \theta\|_{\infty} \geq \varepsilon \right)
    \lesssim_{\gamma} \left( {\color{darkblue}\frac{\sqrt{n}}{t}+\frac{1}{\sqrt{n}}}
    +\frac{\gamma^{\Delta}}{\Delta}+\sqrt{\frac{\Delta}{t}}\right) \varepsilon^{-1}.$
"""

# ╔═╡ 2901a2cb-a2c8-48d9-873e-1a6213d29cad
md"""FIXME
- some kind of method of moments
- main stuff is that we get their limits as functions of $\theta$
-  $1/\sqrt{n}$ comes from the limit of $\hat{v}$ as $n\to \infty$
-  $\sqrt{n}/t$ comes from the limit of $\hat{v}$ as $t\to \infty$, due to the sum of squares
"""

# ╔═╡ 00139eee-0abb-4409-8c96-cd44a6395d79
md"## Community detection"

# ╔═╡ 32e22997-52f9-46e5-8bdc-2b535953a3c6
Columns(md"""
$\mathbb{P}_{\psi}\left( X^{n,i}_{t} = 1 | X^n_{t-1} = x \right) = \mu + \sum_{j\in \mathcal{P}_+} \frac{\gamma}{n}  \psi_{ij} x_j + \sum_{j\in \mathcal{P}_-} \frac{\gamma}{n}\psi_{ij} (1-x_j)$

We expect that $\Sigma^n_{ij} = \operatorname{Cov}_\psi(X^{n,i}_t, X^{n,j}_{t-1})$ is related to $\begin{cases} \psi_{ij}\frac{\gamma}{n} & \text{if } \color{red}{j\in \mathcal{P}_+}\\ -\psi_{ij}\frac{\gamma}{n} & \text{if } \color{blue}{j\in \mathcal{P}_-}\end{cases}$. We prove that:
- the coordinates of $\sigma^n = 1^\top \Sigma^n$ concentrate on two values $\color{red}{\sigma_+}$ and $\color{blue}{\sigma_-}$,
-  $\hat{\sigma}^{n,j}_t=\sum_{i=1}^n \left[\frac{1}{t-1}\sum_{s=2}^t X^{n,i}_s X^{n,j}_{s-1}-\frac{N^{n,i}_{t}}{t}\frac{N^{n,j}_{t}}{t}\right] \to \sigma^{n,j}$ uniformly in $j$.
\

**Detection procedure:** $k$-means ($k=2$) applied to $\hat{\sigma}^n_t$ $\rightsquigarrow$ $\color{red}{\hat{\mathcal{P}}^{n,t}_+}$ and $\color{blue}{\hat{\mathcal{P}}^{n,t}_-}$
""",
    App() do
        fig = Figure(size=(400, 500))

        begin # parameters
            r₊ = 0.5
            μ = 0.25
            γ = 0.5
            p = 0.5
            model = MarkovChainModel(μ, 1.0 - γ, p)
            r₋ = 1.0 - r₊
            D = γ * p * (r₊ - r₋)
            m = (μ + γ * p * r₋) / (1 - D)
            σ₊, σ₋ = (γ * p * m * (1 - m)) .* ((γ^2 * p^2 * (r₊ - r₋)) / (1 - D^2) .+ (1.0, -1.0))

            gltop = GridLayout(fig[1, 1], tellwidth=false)
            Label(gltop[1, 1], L"\sigma_+, \ \sigma_-")
            tog_σs = Toggle(gltop[1, 2], active=false)
            σs_visible = tog_σs.active
            Label(gltop[1, 3], "Reveal")
            tog_color = Toggle(gltop[1, 4], active=false)
            color_visible = tog_color.active
            sg = SliderGrid(fig[2, 1],
                (label=L"n", range=50:10:200, startvalue=50),
                (label=L"t", range=(10:10:300) .^ 2, startvalue=100),
                tellwidth=false)
            n, t = [sl.value for sl in sg.sliders]

            perms = Dict([n => randperm(n) for n in 50:10:200])
        end

        begin # observables
            excitatory = @lift MeanFieldGraph.N2excitatory($n, r₊)
            data = @lift rand(model, $excitatory, 90000)
			xticks_obs = @lift $σs_visible ? ([σ₊, σ₋, 0], [L"\sigma_+", L"\sigma_-", "0"]) : ([],[])
        end

        begin # axis and plot
            σ̂ = @lift MeanFieldGraph.covariance_vector($data[1:$t])
            ax = Axis(fig[3, 1], xlabel=L"\hat{\sigma}^n_t", ylabel="neuron",
                yticks=@lift([1, $n]), xticks=xticks_obs,
                limits=@lift (0.2 .* (-1.0, 1.0), (1 - 0.05 * $n, 1.05 * $n)))

            scatter!(ax, @lift([($σ̂[k], perms[$n][k]) for k in 1:(floor(Int, $n * r₊))]), markersize=5, color=@lift $color_visible ? :red : :black)
            scatter!(ax, @lift([($σ̂[k], perms[$n][k]) for k in (floor(Int, $n * r₊)+1):$n]), markersize=5, color=@lift $color_visible ? :blue : :black)

            vlines!(ax, [σ₊, σ₋], color=[:red, :blue], visible=σs_visible)
        end

        begin # Adjust layout

        end

        DOM.div(fig; style="text-align: center")
    end,
    widths=[65, 35])

# ╔═╡ d3f7f844-14fd-407d-8fd8-94f03afffc2b
md"""
!!! danger "Theorem - CO (2025)"
	If $n\to \infty$ and $\frac{n}{\sqrt{t}}{\small\ln (nt)} \to 0$, then $\mathbb{P}\left(\left\{\hat{{\cal P}}_{+}^{n,t}={\cal P_{+}}\right\}\cap \left\{\hat{{\cal P}}_{-}^{n,t}={\cal P_{-}}\right\} \right) \to 1$.
"""

# ╔═╡ f5f6c19c-b4d1-40b5-8598-fcc66a90dccc
md"""FIXME: fun fact: there are coupling arguments in those papers."""

# ╔═╡ 8d4f6147-758f-4056-9ce2-3baa1513ee9f
md"## Probability of exact recovery and misclassification proportion"

# ╔═╡ 3188c2b6-f08d-4df2-b2cc-8613cb30145b
md"""FIXME
- explain top-left plot
- explain bottom-right plot :
  - difficult when $p$ small
  - exact recovery is way harder than misclassification proportion (especially when $n$ is large)
"""

# ╔═╡ 0d09518c-8e18-45b2-883d-2324c93e68f2
md"### Plot parameters"

# ╔═╡ 42248a46-895a-48e9-8fae-737eb8b63d1d
DOM.div(md"Ribbons' width is two standard deviations.",
    md"""
  -  1000 simulations
  """,
    md"""
    -  $n=50$
    """,
    md"""
    -  $\mu = 0.25$
    """,
    md"""
    -  $\gamma = 0.5$
    """,
    md"""
    -  $p=0.5$
    """,
    md"""
    -  $r_+ = 0.5$
    """)

# ╔═╡ 7b101707-4364-4328-af54-3cea58698ab4
md"## Outlooks"

# ╔═╡ 0c1d8ea0-c116-4392-a0ba-7c11a563a036
md"""FIXME
- overview of my projects organized in the same boxes as before
- Proba : Study strong approximations
  - SPDE
- Julia and the PP.jl
- project with Sophie and Elodie
- statistical models with unknown classes where you provide estimation and classification at the same time
  - types of neurons
  - types of users in social networks
"""

# ╔═╡ 200ae1c8-1521-4f35-98b3-6981d8b7a0d1
html"</br>"

# ╔═╡ 941709a1-2c58-4d29-b128-84212cee668b
md"## Thank you"

# ╔═╡ e65af828-655c-4d4e-a040-4599ee179e89
md"# Appendix"

# ╔═╡ eb5975a0-7732-4033-bc66-e7e429611baf
DOM.div(
    md"- [Age dependent Hawkes process and Refractory density equation](#Age-dependent-Hawkes-process-/-Refractory-Density-Equation)",
    md"- [Sensitivity analysis](#Sensitivity-analysis-of-RDE)",
    md"- [Formal definition of MC-HP](#Multiclass-Hawkes-process)",
    md"- [Convergence of the graph density estimator](#Convergence-of-\hat{\theta})",
    md"- [Optimality of the graph density estimator](#Estimation-optimality)",
    md"- [Optimality of the community detection](#Detection-optimality)",
    md"- [Future of PointProcesses.jl](#PointProcesses.jl)"
)

# ╔═╡ 5f529f1a-cc46-4662-9fbf-0881151e557f
md"## Age dependent Hawkes process / Refractory Density Equation"

# ╔═╡ a2b2253b-89c1-4bc5-a816-ace56a5099f4
html"""<div style="height: 100px;"></div>"""

# ╔═╡ 301274ed-109a-44d7-a869-97a49f2bb30b
md"### Plot parameters"

# ╔═╡ c7b7ff57-5671-447d-9f15-78236f9ffa1a
DOM.div(md"""
-  $n=100$
""",
    md"""
    -  $f(a, \xi) = \phi(\xi) \mathbf{1}_{[a^*, \infty)}(a)$
    """,
    md"""
    -  $\phi(\xi) = 1/2 + 10\xi^2 / (\xi^2 + 1)$
    """,
    md"""
    -  $h(t) = \alpha e^{-\alpha t}$ with $\alpha = 10$
    """,
    md"""
    -  $u^{\rm in}(a) = e^{-(a -1)}\mathbf{1}_{[1,\infty)}(a)$
    """)

# ╔═╡ 6273cbc1-ca15-43ad-bd66-a7d89232201b
md"## Sensitivity analysis of RDE"

# ╔═╡ b3275370-4e8c-41fe-bd70-f6bc7eb04a5b
Columns(DOM.div(md"### Bio realistic parameters",
        md"""
  -  $a^* = 5$ms corresponds to $1/a^* = 200$Hz
  """,
        md"""
      -  $\mu = 2$Hz
      """,
        md"""
      -  $\mu a^* = 0.01$
      """),
    DOM.div(PlutoUI.LocalResource("assets/RDE-sensitivity_max.png", :width => 420); style="text-align: center"))

# ╔═╡ be991ec2-e0e5-40f5-8451-c8a1a850e605
md"## Multiclass Hawkes process"

# ╔═╡ c019a231-464e-4a2c-95cf-f806cbea5c4a
md"## Convergence of $\hat{\theta}$"

# ╔═╡ 1eca6154-14ce-4a73-8983-e686997b1b55
html"""<div style="height: 100px;"></div>"""

# ╔═╡ 7ecbde43-ef95-4ed9-942d-149a346ca781
md"### Plot parameters"

# ╔═╡ 5a917ed1-cfc1-40f7-838e-809bdd5a4db3
DOM.div(md"""
-  100 simulations
""",
    md"""
    -  $n=500$
    """,
    md"""
    -  $\mu = 0.25$
    """,
    md"""
    -  $\gamma = 0.5$
    """,
    md"""
    -  $p=0.5$
    """,
    md"""
    -  $r_+ = 0.5$
    """)

# ╔═╡ 09b7fe1f-4994-42a9-9db0-5f2472262124
md"## Estimation optimality"

# ╔═╡ 8062da1a-5b93-4510-918b-dfc8ca767b93
md"[Link to the upper-bound](FIXME)"

# ╔═╡ 28ac8b5f-06c3-4923-b142-24fd7cf8685b
md"""
-  $\Psi \sim \operatorname{Bin}(n,p)$ with unknown $p$. 
- Let $K$ be a known value and define $\gamma(p)= K/p$ and $m=1/2+K$
- Simplified version of the model: $\mu = 1/2$ is known, $\mathcal{P}_- = \emptyset$, no temporal or spatial dependence.
"""

# ╔═╡ 932dfc39-958a-4df7-acd0-bcebe912b48e
md"""
### Setting 1
- Given $\Psi = \psi$, we observe a $n$-sample of r.v. $B \sim \operatorname{Bin}(t,1/2+\gamma(p)\psi/n)$.

!!! danger ""
	The estimator of $1/p$ from the method of moments is unbiased and its standard deviation is of order $\frac{\sqrt{n}}{t}+\frac{1}{\sqrt{n}}$.
"""

# ╔═╡ 37c5f7e7-5307-4493-8c58-7fce21c969bd
md"""
### Setting 2
- Given $\Psi = \psi$, we observe a $n\lfloor t/2 \rfloor$-sample of r.v. $B \sim \operatorname{Bin}(2,1/2+\gamma(p)\psi/n)$.

!!! danger ""
	The standard deviation of any unbiased estimator of $1/p$ is of order at least $\frac{\sqrt{n}}{\sqrt{t}}$.
"""

# ╔═╡ 0b48d6b6-e4a0-4355-9722-ce4ccd30a560
md"## Detection optimality"

# ╔═╡ 8799c616-5967-49c8-b62a-1aef442bb187
@markdown(
    """
    - Estimation/classification in turn: fix ``r_+`` ``\\to`` estimate ``\\theta`` ``\\to`` plug into ``(\\sigma_+, \\sigma_-)`` ``\\to`` estimate ``r_+`` ``\\to`` ``\\dots`` ``\\quad`` <span id="small-gray" style="float:right">(thanks Jeff)</span>
    """
)

# ╔═╡ 41b02bcf-40cf-4234-9bdc-498a29ae6206
DOM.div(md"""
!!! tip "Lower bound (application of Assouad's Lemma)"
	If $\frac{t}{n} \to 0$ then the mean misclassification proportion $\not \rightarrow 0$.
""",
    md"""### Spectral procedure

  -  $(\Sigma^n)^\top \Sigma^n$ is close to a rank one matrix
  - Its eigenvector is $\sigma_+ \mathbf{1}_{\mathcal{P}_+} + \sigma_- \mathbf{1}_{\mathcal{P}_-}$
  """,
    md"""
    - Clustering based on the largest eigenvector of $(\hat{\Sigma}^n)^\top \hat{\Sigma}^n$
    """,
    md"""
!!! tip "Upper-bound"
	The mean misclassification proportion of the procedure is of the order $\frac{n}{t}$ (up to $\log$ terms) with high probability.

- Concentration inequalities in terms of spectral norm
- Coupling similar to Ost & Reynaud-Bouret (2020)
- Pertubation theory for matrices
""")

# ╔═╡ 355cacab-cbed-49a2-8776-9ea4b341c01d
md"## PointProcesses.jl"

# ╔═╡ ea823dca-9737-4237-911b-d8ba9debf830
@markdown("""
- Benchmark simulation algorithms of linear Hawkes processes
</br>
		  
- Implement KMT coupling in both directions
</br>
	  
- Study the computational tradeoff involved in thinning simulation
</br>
		  
- Use binary heaps for sparse Hawkes processes
</br>
		  
- Conformal prediction for point processes
""")

# ╔═╡ c71a8f79-762a-4568-93f0-ea1a186eaec8
md"# Pluto Appendix"

# ╔═╡ dd3c0f63-8bf7-41d2-a958-a681894640ea
slideresolution = (0.96, 0.81) .* monitorresolution;

# ╔═╡ 7a9f4d5f-76a4-4484-b0b4-661ceed69163
slideresolution

# ╔═╡ 84e32427-cc6b-4b35-8a67-a243d8ec8add
App() do
    # Constants for layout and geometry
    FIGURE_RESOLUTION = 1.04 .* slideresolution
    TRIANGLE_RADIUS = 2.0             # Distance from the middle box to the outer boxes
    BOX_HEIGHT = 1.0                  # Height of box
    BOX_WIDTH = 16 / 9 * BOX_HEIGHT   # Width of box
    CORNER_RADIUS = 0.1               # Radius of box corners
    ZOOM_PADDING = 1.1                # Padding factor around zoomed box
    ANIMATION_DURATION = 2.0          # Animation duration in seconds
    FRAMES_PER_SECOND = 60            # Animation frames per second
    EASING_FUNCTION = t -> 3t^2 - 2t^3  # Smooth step easing (cubic Hermite)
    BOX_LABELS = ("""
        Point processes
        Neurosciences
        """,
        "Probability",
        "PDE",
        "Statistics")

    # Function to draw a rounded rectangle box
    function RoundedRectangle(center::Point2, width::Real, height::Real, radius::Real)
        # Clamp radius to valid range
        r = clamp(Float64(radius), 0.0, min(width, height) / 2)
        # retrieve coordinates
        xc, yc = center
        # compute bottomleft coordinates
        x, y = (xc - width / 2, yc - height / 2)

        # Zero radius? Use a standard rectangle
        if r < 1e-10
            return Rect(x, y, width, height)
        end

        # Bezier control point offset for quarter circle approximation
        # This magic number is 4/3 * (sqrt(2) - 1) ≈ 0.5522847498
        k = 0.5522847498 * r

        # Build the path (clockwise from top edge)
        return BezierPath([
            MoveTo(Point2f(x + width - r, y)),
            LineTo(Point2f(x + r, y)),
            CurveTo(Point2f(x + r - k, y), Point2f(x, y + r - k), Point2f(x, y + r)),
            LineTo(Point2f(x, y + height - r)),
            CurveTo(Point2f(x, y + height - r + k), Point2f(x + r - k, y + height), Point2f(x + r, y + height)),
            LineTo(Point2f(x + width - r, y + height)),
            CurveTo(Point2f(x + width - r + k, y + height), Point2f(x + width, y + height - r + k), Point2f(x + width, y + height - r)),
            LineTo(Point2f(x + width, y + r)),
            CurveTo(Point2f(x + width, y + r - k), Point2f(x + width - r + k, y), Point2f(x + width - r, y)),
            ClosePath()
        ])
    end

    # Animation function with lock to prevent concurrent animations
    function animate_limits!(target_x, target_y; zoom_in=true)
        global is_animating, current_xlims, current_ylims

        if is_animating
            return
        end

        # Does nothing if already at target
        if target_x == current_xlims && target_y == current_ylims
            return
        end

        is_animating = true

        # Calculate animation parameters
        total_frames = ceil(Int, ANIMATION_DURATION * FRAMES_PER_SECOND)
        middle_frame = div(total_frames, 2)
        start_xlims = current_xlims
        start_ylims = current_ylims

        for frame in 1:total_frames
            # Check if figure is still open
            if !isopen(fig.scene)
                is_animating = false
                return
            end

            t_raw = frame / total_frames
            t_switch_in = 0.6
            t_switch_out = 1.0 - t_switch_in

            if zoom_in
                alpha_background[] = EASING_FUNCTION(t_raw / t_switch_in) * (t_raw <= t_switch_in) + EASING_FUNCTION(1.0 - (t_raw - t_switch_in) / (1.0 - t_switch_in)) * (t_raw > t_switch_in)
                notify(alpha_background)
                t_raw == t_switch_in && (visibles_labels[] = false; notify(visibles_labels))
                t_raw == t_switch_in && (visible_images[] = true; notify(visible_images))
            end
            if !zoom_in
                alpha_background[] = EASING_FUNCTION(t_raw / t_switch_out) * (t_raw <= t_switch_out) + EASING_FUNCTION(1.0 - (t_raw - t_switch_out) / (1.0 - t_switch_out)) * (t_raw > t_switch_out)
                notify(alpha_background)
                t_raw == t_switch_out && (visibles_labels[] = true; notify(visibles_labels))
                t_raw == t_switch_out && (visible_images[] = false; notify(visible_images))
            end

            # Apply easing function for smooth acceleration/deceleration
            t_smooth = EASING_FUNCTION(t_raw)

            # Interpolate limits
            new_xlims = (
                start_xlims[1] + (target_x[1] - start_xlims[1]) * t_smooth,
                start_xlims[2] + (target_x[2] - start_xlims[2]) * t_smooth
            )
            new_ylims = (
                start_ylims[1] + (target_y[1] - start_ylims[1]) * t_smooth,
                start_ylims[2] + (target_y[2] - start_ylims[2]) * t_smooth
            )

            # Update axis limits
            xlims!(ax, new_xlims)
            ylims!(ax, new_ylims)

            # Update current state
            current_xlims = new_xlims
            current_ylims = new_ylims

            # Frame timing
            sleep(1 / FRAMES_PER_SECOND)
        end

        # Ensure final position is exact
        xlims!(ax, target_x)
        ylims!(ax, target_y)
        current_xlims = target_x
        current_ylims = target_y

        is_animating = false
    end

    function animate_limits_2steps!(target_x, target_y)
        global is_animating, current_xlims, current_ylims

        # Does nothing if already at target
        if target_x == current_xlims && target_y == current_ylims
            return
        end

        @async begin
            # Reset view
            animate_limits!(initial_xlims, initial_ylims; zoom_in=false)
            # Move view
            animate_limits!(target_x, target_y; zoom_in=true)
        end
    end

    begin # Initialize figure
        # Create figure and axis
        fig = Figure(size=FIGURE_RESOLUTION)
        ax = Axis(fig[1, 2], aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)

        # Centers of circles
        centers = [
            Point2f(0, 0),
            Point2f(0, TRIANGLE_RADIUS),
            Point2f(-TRIANGLE_RADIUS, -0.7 * TRIANGLE_RADIUS),
            Point2f(TRIANGLE_RADIUS, -0.7 * TRIANGLE_RADIUS)
        ]

        # Strokecolors and labels for the circles
        colors = [:blue, :green, :red, :purple]
        visibles_labels = Observable(true)
        visible_images = Observable(false)
        alpha_background = Observable(0.0)

        for i in 1:4
            # Box content as image
            xpad = BOX_WIDTH / 2
            ypad = BOX_HEIGHT / 2
            img_x = centers[i][1] .+ (-1.0, +1.0) .* xpad
            img_y = centers[i][2] .+ (-1.0, +1.0) .* ypad
            img = load("assets/box" * string(i) * ".png")
            image!(ax, img_x, img_y, rotr90(img), visible=visible_images)

            # Box label
            text!(ax, centers[i],
                text=BOX_LABELS[i],
                align=(:center, :center),
                fontsize=28,
                color=:black,
                font=:bold,
                visible=visibles_labels)

            # Box background to deal with transparency
            poly!(ax, RoundedRectangle(centers[i], BOX_WIDTH, BOX_HEIGHT, CORNER_RADIUS),
                color=:white,
                alpha=alpha_background)

            # Box contour
            poly!(ax, RoundedRectangle(centers[i], BOX_WIDTH, BOX_HEIGHT, CORNER_RADIUS),
                color=:transparent,
                strokecolor=Makie.Colors.JULIA_LOGO_COLORS[colors[i]],
                strokewidth=4)
        end

        # Set initial limits to show all boxes
        x_extent = TRIANGLE_RADIUS + BOX_WIDTH / 2
        y_extent = TRIANGLE_RADIUS + BOX_HEIGHT / 2
        initial_xlims = (-x_extent * 1.1, x_extent * 1.1)
        initial_ylims = (-y_extent * 1.1, y_extent * 1.1)
        xlims!(ax, initial_xlims)
        ylims!(ax, initial_ylims)

        # Create selection buttons
        buttongrid = fig[1, 1] = GridLayout(tellheight=false)
        buttonlabels = ["Show all", "PP/Neuro", "Probability", "PDE", "Statistics"]
        buttons = buttongrid[1:5, 1] = [Makie.Button(fig, label=buttonlabels[i], width=Fixed(100)) for i in 1:5]

        # Animation state
        global is_animating = false
        global current_xlims = initial_xlims
        global current_ylims = initial_ylims
    end

    # Connect buttons to zoom action
    for buttonid in 1:5
        on(buttons[buttonid].clicks) do _
            if buttonid == 1
                target_x = initial_xlims
                target_y = initial_ylims
            else
                circle_num = buttonid - 1
                center = centers[circle_num]

                xpad = BOX_WIDTH * ZOOM_PADDING / 2
                ypad = BOX_HEIGHT * ZOOM_PADDING / 2
                target_x = (center[1] - xpad, center[1] + xpad)
                target_y = (center[2] - ypad, center[2] + ypad)
            end

            # Start animation
            animate_limits_2steps!(target_x, target_y)
        end
    end

    # Adjust layout for better appearance
    colsize!(fig.layout, 1, Fixed(150))

    DOM.div(fig; style="text-align: center")
end

# ╔═╡ 78de54a6-6f9a-4a5c-8121-48c9c2a80020
DOM.div(PlutoUI.LocalResource("assets/voltage-trace.png", :width => slideresolution[1] / 2), style="text-align:center")

# ╔═╡ 655cddf9-34f5-41c5-88b8-d095c3dbce5a
App() do
    function continuous2discrete(h::History, nframes::Int)
        frame = 1
        tmin = min_time(h)
        tmax = max_time(h)
        discrete_grid = range(tmin, tmax, nframes + 1)[2:end]
        discrete_t = discrete_grid[1]
        times = event_times(h)
        marks = event_marks(h)
        Nneur = max_mark(h)

        output = zeros(Bool, Nneur, nframes)
        for (i, t) in enumerate(times)
            while t > discrete_t
                frame += 1
                discrete_t = discrete_grid[frame]
            end

            output[marks[i], frame] = 1
        end

        return output
    end

    function fadedots(; fps::Int=60, fade_rate=0.01)
        while any(alphas[] .!= 0)
            alphas[] .= clamp.((alphas[] .- fade_rate), 0.0, 1.0)
            sleep(1 / fps)
        end
    end

    function spikesanimation(spikes::Matrix; fps::Int=60, fade_rate=0.01)
        # Reset scatter plot
        stopflag[] = true
        fadedots(fps=fps, fade_rate=fade_rate) # fade existing dots
        sleep(0.1)

        nframes = size(spikes)[2]
        # Start new animation
        stopflag[] = false
        @async begin
            for frame in 1:nframes
                stopflag[] == true && break

                # Fade existing dots
                alphas[] .= clamp.((alphas[] .- fade_rate), 0.0, 1.0)

                # Add new dots
                alphas[][findall(spikes[:, frame] .== true)] .= 1.0

                notify(alphas)
                sleep(1 / fps)
            end
            # Fade remaining dots
            for _ in 1:ceil(Int, 1 / fade_rate)
                alphas[] .= clamp.((alphas[] .- fade_rate), 0.0, 1f0)
                notify(alphas)
                sleep(1 / fps)
            end
        end
    end

    function spikesanimation(; fps::Int=60, fade_rate=0.01)
        stopflag[] = true
        fadedots(fps=fps, fade_rate=fade_rate) # fade existing dots
        sleep(0.1)

        # Start new animation
        stopflag[] = false
        @async begin
            while !stopflag[]
                # Fade existing dots
                alphas[] .= clamp.(alphas[] .- fade_rate, 0.0, 1.0)

                if rand() < 0.05
                    alphas[][rand(DiscreteUniform(1, Nneur))] = 1.0
                end
                notify(alphas)
                sleep(1 / fps)
            end
            # Fade remaining dots
            for _ in 1:ceil(Int, 1 / fade_rate)
                alphas[] .= clamp.((alphas[] .- fade_rate), 0.0, 1f0)
                notify(alphas)
                sleep(1 / fps)
            end
        end
    end

    begin # Load simulated data 
        # Load positions of neurons corresponding to a SHP network with 900 neurons
        positions = load("../data/SHP-animation.jld2", "positions")
        Nneur = length(positions)

        # Load spikes
        ## Initialize a vector of spike data
        spikes = Matrix{Bool}[]
        ## Load spike history made from a ADHP network with 900 neurons
        ### generated by animation-ADHP.jl
        h_adhp = load("../data/ADHP-animation.jld2", "history")
        push!(spikes, continuous2discrete(h_adhp, 2000))
        ## Load spike history made from a SHP network with 900 neurons
        ### generated by animation-SHP.jl
        h_shp = load("../data/SHP-animation.jld2", "history")
        push!(spikes, continuous2discrete(h_shp, 2000))
        ## Load spike history made from a 2CHP network with 900 neurons
        ### generated by animation-2CHP.jl
        h_2chp = load("../data/2CHP-animation.jld2", "history")
        push!(spikes, continuous2discrete(h_2chp, 2000))
        ## Load discrete spikes made from a MCRE network with 900 neurons
        ### generated by animation-MCRE.jl
        push!(spikes, load("../data/MCRE-animation.jld2", "spikes"))
    end

    begin
        buttonlabels = ["Age dependent Hawkes process",
            "Spatial Hawkes process",
            "Multi class Hawkes process",
            "Binary Markov chain"]
    end

    begin # Animation setup
        # Figure setup
        fig = Figure(size=slideresolution, backgroundcolor=:black)
        buttongrid = fig[1, 1] = GridLayout(tellheight=false)
        startbutton = buttongrid[1, 1] = Makie.Button(fig, label="Start", fontsize=30)

        # Observables
        alphas = Observable(zeros(Nneur))
        stopflag = Observable(false)

        on(startbutton.clicks) do _
            # Figure update
            ## Remove start button
            empty!(startbutton.blockscene)
            delete!(startbutton)
            ## Initialize scatter plot
            ax = Axis(fig[1, 2], backgroundcolor=:black, aspect=AxisAspect(1), topspinecolor=:white, leftspinecolor=:white, bottomspinecolor=:white, rightspinecolor=:white)
            scatter!(ax, positions, color=@lift(RGBAf.(1, 1, 1, $alphas)), markersize=8)
            ## Initialize buttons
            buttons = buttongrid[1:4, 1] = [Makie.Button(fig, label=buttonlabels[i], fontsize=20) for i in 1:4]

            fps = 60
            fade_rate = 0.01
            spikesanimation(; fps=fps, fade_rate=fade_rate)

            for buttonid in 1:4
                on(buttons[buttonid].clicks) do _
                    spikesanimation(spikes[buttonid]; fps=fps, fade_rate=fade_rate)
                end
            end
        end
    end

    DOM.div(fig; style="text-align: center")
end

# ╔═╡ c49ac056-ab65-44ea-b37c-2e9d6426399e
Columns(DOM.div(md"""
- Intensity of ``N^i`` (w.r.t. filtration ``\mathcal{F}_t``):
$\lambda^i_t = \lim_{\mathrm{d}t \to 0} \frac{1}{\mathrm{d}t} \mathbb{P}\left( N^i \cap (t, t + \mathrm{d}t] =1 \middle| \mathcal{F}_t \right).$
""",
        html"</br>",
        @markdown("""
        - Models <span id="red">excitatory</span>/<span id="blue">inhibitory</span> behaviors.
        """)),
    App() do
        fig = Figure(size=(0.4 * slideresolution[1], 0.35 * slideresolution[2]))
        N1 = [0.43, 2.54, 5.43, 6.12, 8.34]
        y1 = 6
        N2 = [0.54, 0.76, 2.84, 2.96, 3.12, 5.65, 5.90, 6.15, 6.45, 6.87, 8.54, 8.90]
        y2 = 5
        N3 = [1.32, 3.45, 6.25, 9.40]
        y3 = 2
        N4 = [0.5, 1.03, 2.55, 3.32, 4.65, 5.12, 5.34, 5.98, 7.25, 7.56, 8.01, 8.45, 8.61, 8.96]
        y4 = 1

        axtop = Axis(fig[1, 1], yticks=([y1, y2, y3, y4], [L"i", L"j", L"i", L"j"]), xticks=([0], [""]), yticklabelsize=20)
        hidedecorations!(axtop, ticklabels=false)
        hidespines!(axtop)
        scatter!(axtop, N1, fill(y1, length(N1)), color=:black)
        scatter!(axtop, N2, fill(y2, length(N2)), color=:black)
        scatter!(axtop, N3, fill(y3, length(N3)), color=:black)
        scatter!(axtop, N4, fill(y4, length(N4)), color=:black)
        for x in N1
            arrows2d!(axtop, Point(x, y1), Vec(0.4, -1), lengthscale=0.6, align=-0.2, color=:red)
        end
        for x in N3
            arrows2d!(axtop, Point(x, y3), Vec(0.4, -1), lengthscale=0.6, align=-0.2, color=:blue)
        end

        DOM.div(fig)
    end,
    widths=[60, 40])

# ╔═╡ 357a9f03-912a-4ea0-9d0c-39ff734ebf01
App() do
    fig = Figure(size=(slideresolution[1], slideresolution[2] - 300))

    # Compute a 2d Hawkes process and its intensities via thinning
    function thinning(h_poisson, μ, α, β, dt)
        tmin = min_time(h_poisson)
        tmax = max_time(h_poisson)
        times_poisson = event_times(h_poisson)
        marks_poisson = event_marks(h_poisson)

        h = History(; times=Float64[], marks=Int[], tmin=tmin, tmax=tmax)

        λ = μ
        t = tmin
        Ξ = fill(0, length(μ))
        λs = λ
        ts = [t]
        i = 1
        while t < tmax
            τ = (i < length(times)) ? times_poisson[i] - t : Inf
            if τ > dt
                Ξ *= exp(-α * dt)
                λ = max.(0, μ .+ Ξ)
                λs = hcat(λs, λ)
                t += dt
                push!(ts, t)
            else
                Ξ *= exp(-α * τ)
                λ = max.(0, μ .+ Ξ)
                λs = hcat(λs, λ)
                t += τ
                push!(ts, t)

                j, u = marks[i]
                m = Int(j)
                if u < λ[m]
                    push!(h, t, m)
                    Ξ += β[:, m]
                    λ = max.(0, μ .+ Ξ)

                    λs = hcat(λs, λ)
                    push!(ts, t)
                end

                i += 1
            end
        end
        return ts, λs, h
    end

    intensity_bound = 5.0
    begin # time grid
        tmin = 0.0
        tmax = 5.0
        dt = 0.05
    end

    begin # simulation of two 2d Poisson for thinning
        seed!(2)
        N = rand(Poisson(float(2 * intensity_bound * (tmax - tmin))))
        times = sort!(rand(Uniform(tmin, tmax), N))
        marks = [rand(product_distribution(DiscreteUniform(1, 2), Uniform(0.0, intensity_bound))) for n in 1:N]
        h_poisson = History(; times=times, marks=marks, tmin=tmin, tmax=tmax)
    end

    begin # parameters 
        sg = SliderGrid(fig[1, 1],
            (label=L"\mu_1", range=0.8:0.01:0.95, startvalue=0.9),
            (label=L"\mu_2", range=1.2:0.01:1.4, startvalue=1.3),
            (label=L"\alpha", range=1.0:0.01:5.0, startvalue=4.0),
            (label=L"w_{11}", range=-1.0:0.01:1.0, startvalue=0.5),
            (label=L"w_{12}", range=-1.0:0.01:1.0, startvalue=-0.2),
            (label=L"w_{21}", range=-1.0:0.01:1.0, startvalue=-0.5),
            (label=L"w_{22}", range=-1.0:0.01:1.0, startvalue=0.3),
            tellheight=false)
        obs_parameters = sg.sliders
        μ = @lift [$(obs_parameters[1].value), $(obs_parameters[2].value)]
        α = obs_parameters[3].value
        β = @lift [$(obs_parameters[4].value) $(obs_parameters[6].value);
            $(obs_parameters[5].value) $(obs_parameters[7].value)]
    end

    begin # Observables
        obs_thinned = @lift thinning(h_poisson, $μ, $α, $β, dt)
        ts = @lift $obs_thinned[1]
        λs = @lift $obs_thinned[2]
        h = @lift $obs_thinned[3]

        spikes1 = @lift event_times($h)[event_marks($h).==1]
        spikes2 = @lift event_times($h)[event_marks($h).==2]
    end

    begin # Plot
        ax = Axis(fig[1, 2], ylabel=L"\lambda", xlabel=L"t", title="Bivariate Hawkes process")

        lines!(ax, ts, (@lift $λs[1, :]), color=:darkblue, label=L"\lambda^1_t")
        lines!(ax, ts, (@lift $λs[2, :]), color=:orange, label=L"\lambda^2_t")
        hlines!(ax, (@lift $μ[1]), color=:darkblue, linestyle=:dot)
        hlines!(ax, (@lift $μ[2]), color=:orange, linestyle=:dot)
        scatter!(ax, spikes1, (@lift 0.0 * $spikes1), color=:darkblue, markersize=15, label=L"N^1")
        scatter!(ax, spikes2, (@lift 0.0 * $spikes2), color=:orange, markersize=15, label=L"N^2")

        axislegend(ax)
    end

    begin # Adjust layout
        colsize!(fig.layout, 1, Fixed(300))
    end

    DOM.div(fig; style="text-align: center")
end

# ╔═╡ 0915e18b-0f4b-4fdd-8ab2-ceb26a0eff7b
Columns(md"""
!!! danger "Thinning"
	Let $\Pi$ be a 2d-Poisson process with intensity $1$. Then, $(\lambda_t)_{t\geq 0}$ is the intensity of $N$ defined by
	
	$N(\mathrm{d}t) = \int_{0}^{\infty} \mathbf{1}_{[0,\lambda_t]}(z) \, \Pi(\mathrm{d} t, \mathrm{d} z).$
""", App() do
        fig = Figure(size=(700, slideresolution[2] / 2))

        intensity_bound = 1.5
        begin # time grid
            tmin = 0.0
            tmax = 5.0
        end

        begin # user-defined Poisson for thinning
            h_poisson = History(;
                times=[0.75, 1.32, 2.01, 2.87, 3.78, 4.54],
                marks=[1.23, 0.95, 0.43, 0.67, 1.18, 0.87],
                tmin=tmin, tmax=tmax
            )
            spikes12 = [2.01, 2.87]
            marks12 = [0.43, 0.67]
            spikes1 = [1.32]
            marks1 = [0.95]
            spikes2 = [4.54]
            marks2 = [0.87]

            # Defines points and vectors to draw arrows
            crosspad = 0.05
            crosses12 = Point2f.(spikes12, marks12 .- crosspad)
            vec12 = Vec2f.(0, -(marks12 .- 2 * crosspad))
            crosses1 = Point2f.(spikes1, marks1 .- crosspad)
            vec1 = Vec2f.(0, -(marks1 .- 2 * crosspad))
            crosses2 = Point2f.(spikes2, marks2 .- crosspad)
            vec2 = Vec2f.(0, -(marks2 .- 2 * crosspad))
        end

        begin # parameters 
            sg = SliderGrid(fig[2, 1],
                (label="t", range=tmin:0.01:tmax, startvalue=tmin), tellwidth=false)
            t_animation = sg.sliders[1].value
        end

        begin # observables
            x_range_1 = @lift [tmin, min($t_animation, tmin + (tmax - tmin) / 2)]
            x_range_2 = @lift [min($t_animation, tmin + (tmax - tmin) / 2), min($t_animation, tmax)]

            spikes12_animated = @lift spikes12[spikes12.<$t_animation]
            spikes1_animated = @lift spikes1[spikes1.<$t_animation]
            spikes2_animated = @lift spikes2[spikes2.<$t_animation]
        end

        begin # Plot layout parameters
            alpha_intensity = 0.05
            marksize = 13
            ypad = 0.1
        end

        begin # Plot
            ax = Axis(fig[1, 1], ylabel=L"\lambda", xlabel="t", title="SDE coupling", xticks=[0, 2.5, 5])
            limits!(ax, (nothing, nothing), (0.0 - ypad, intensity_bound + ypad))

            # vertical line to represent t_animation
            vlines!(ax, t_animation, color=:red, alpha=0.5)

            # transparent rectangles for area under the intensity
            poly!(ax, Rect(tmin, 0.0, tmin + (tmax - tmin) / 2, 1.0), color=:darkblue, alpha=alpha_intensity)
            poly!(ax, Rect(tmin + (tmax - tmin) / 2, 0.0, tmin + (tmax - tmin) / 2, 0.8), color=:darkblue, alpha=alpha_intensity)
            poly!(ax, Rect(tmin, 0.0, tmin + (tmax - tmin) / 2, 0.8), color=:orange, alpha=alpha_intensity)
            poly!(ax, Rect(tmin + (tmax - tmin) / 2, 0.0, tmin + (tmax - tmin) / 2, 1.0), color=:orange, alpha=alpha_intensity)

            # transparent lines for future intensity
            lines!(ax, [tmin, tmin + (tmax - tmin) / 2], [1.0, 1.0], color=:darkblue, alpha=8 * alpha_intensity)
            lines!(ax, [tmin + (tmax - tmin) / 2, tmax], [0.8, 0.8], color=:darkblue, alpha=8 * alpha_intensity)
            lines!(ax, [tmin, tmin + (tmax - tmin) / 2], [0.8, 0.8], color=:orange, alpha=7 * alpha_intensity)
            lines!(ax, [tmin + (tmax - tmin) / 2, tmax], [1.0, 1.0], color=:orange, alpha=7 * alpha_intensity)

            # crosses for the 2d Poisson
            scatter!(ax, h_poisson.times, h_poisson.marks, color=:black, marker=:xcross, markersize=marksize, label=L"\Pi")

            # solid lines for past intensity
            lines!(ax, x_range_1, [1.0, 1.0], color=:darkblue, label=L"\lambda^1(t)")
            lines!(ax, x_range_2, [0.8, 0.8], color=:darkblue)
            lines!(ax, x_range_1, [0.8, 0.8], color=:orange, label=L"\lambda^2(t)")
            lines!(ax, x_range_2, [1.0, 1.0], color=:orange)

            # dots for the projection of thinned points
            scatter!(ax, spikes12_animated, (@lift 0.0 * $spikes12_animated), color=:brown, markersize=marksize)
            scatter!(ax, spikes1_animated, (@lift 0.0 * $spikes1_animated), color=:darkblue, markersize=marksize, label=L"N^1")
            scatter!(ax, spikes2_animated, (@lift 0.0 * $spikes2_animated), color=:orange, markersize=marksize, label=L"N^2")

            # projection arrows
            arrows2d!(ax, crosses12[1], vec12[1], minshaftlength=0, color=:brown, visible=@lift spikes12[1] < $t_animation)
            arrows2d!(ax, crosses12[2], vec12[2], minshaftlength=0, color=:brown, visible=@lift spikes12[2] < $t_animation)
            arrows2d!(ax, crosses1, vec1, minshaftlength=0, color=:darkblue, visible=@lift spikes1[1] < $t_animation)
            arrows2d!(ax, crosses2, vec2, minshaftlength=0, color=:orange, visible=@lift spikes2[1] < $t_animation)

            axislegend(ax, position=:ct, orientation=:horizontal)
        end

        DOM.div(fig; style="text-align: center")
    end
)

# ╔═╡ 8d01edc8-71ee-4431-a68b-3a3745428176
Columns(md"""
!!! danger "Time change"
	If $N$ has intensity $(\lambda_t)_{t\geq 0}$, then the time changed process $\Pi_\tau = N_{\Lambda^{-1}_{\tau}}$ is a 1d-Poisson process with intensity 1,
	where $\Lambda^{-1}_{\tau} := \inf\{t\geq 0, \Lambda_t \geq \tau\}$ and $\Lambda_t := \int_{0}^{t} \lambda_s \mathrm{d} s$.
""", App() do
        fig = Figure(size=(700, slideresolution[2] / 2))

        begin # time grid
            tmin::Float64 = 0.0
            tmax::Float64 = 5.0
            ts = tmin:0.01:tmax
        end

        begin # user-defined Poisson for thinning
            h_poisson = History(;
                times=[0.75, 2.01, 2.87, 3.78, 4.94],
                marks=fill(nothing, 5),
                tmin=tmin, tmax=tmax
            )
        end

        begin # parameters 
            sg = SliderGrid(fig[2, 1],
                (label="t", range=tmin:0.01:tmax, startvalue=tmin), tellwidth=false)
            t_animation = sg.sliders[1].value

            Λ1(t) = (t - tmin) * (t <= tmin + (tmax - tmin) / 2) + ((tmax - tmin) / 2 + 0.8 * (t - (tmax - tmin) / 2)) * (tmin + (tmax - tmin) / 2 < t)
            Λ1_inv(τ) = τ * (τ <= (tmax - tmin) / 2) + ((tmax - tmin) / 2 + (τ - (tmax - tmin) / 2) / 0.8) * ((tmax - tmin) / 2 < τ)
            Λ2(t) = 0.8 * (t - tmin) * (0 < t <= tmin + (tmax - tmin) / 2) + (0.8 * (tmax - tmin) / 2 + (t - (tmax - tmin) / 2)) * (tmin + (tmax - tmin) / 2 < t <= tmax)
            Λ2_inv(τ) = τ / 0.8 * (τ <= 0.8 * (tmax - tmin) / 2) + ((tmax - tmin) / 2 + (τ - 0.8 * (tmax - tmin) / 2)) * (0.8 * (tmax - tmin) / 2 < τ)
        end

        begin # observables
            t_range = @lift tmin:0.01:$t_animation

            visibles1 = @lift h_poisson.times .< Λ1($t_animation)
            visibles2 = @lift h_poisson.times .< Λ2($t_animation)

            spikes1_animated = @lift Λ1_inv.(h_poisson.times[$visibles1])
            spikes2_animated = @lift Λ2_inv.(h_poisson.times[$visibles2])
        end

        begin # Plot layout parameters
            alpha_intensity = 0.3
            marksize = 13
            arrowspad = 0.2
        end

        begin # Plot
            ax = Axis(fig[1, 1], ylabel="τ", xlabel="t", title="Time change coupling", xticks=[0, 2.5, 5])

            # vertical line to represent t_animation
            vlines!(ax, t_animation, color=:red, alpha=0.5)

            # transparent lines for future cumulative intensity
            lines!(ax, ts, Λ1.(ts), color=:darkblue, alpha=alpha_intensity)
            lines!(ax, ts, Λ2.(ts), color=:orange, alpha=alpha_intensity)

            # dots for the 1d Poisson on y axis
            scatter!(ax, fill(0.0, length(h_poisson)), h_poisson.times, color=:black, markersize=marksize, label=L"\Pi")

            # solid lines for past cumulative intensity
            lines!(ax, t_range, (@lift Λ1.($t_range)), color=:darkblue, label=L"\Lambda^1(t)")
            lines!(ax, t_range, (@lift Λ2.($t_range)), color=:orange, label=L"\Lambda^2(t)")

            # dots for the time changed points
            scatter!(ax, spikes1_animated, (@lift 0.0 * $spikes1_animated), color=:darkblue, markersize=marksize, label=L"N^1")
            scatter!(ax, spikes2_animated, (@lift 0.0 * $spikes2_animated), color=:orange, markersize=marksize, label=L"N^2")

            # projection arrows
            for i in 1:4
                lines!(ax, [0.0, Λ1_inv(h_poisson.times[i])], fill(h_poisson.times[i], 2), color=:darkblue, alpha=0.5, visible=@lift $visibles1[i])
                arrows2d!(ax, Point2f(Λ1_inv(h_poisson.times[i]), h_poisson.times[i]), Vec2f(0, -h_poisson.times[i] + arrowspad), minshaftlength=0, color=:darkblue, alpha=0.5, visible=@lift $visibles1[i])
                lines!(ax, [0.0, Λ2_inv(h_poisson.times[i])], fill(h_poisson.times[i], 2), color=:orange, alpha=0.5, visible=@lift $visibles2[i])
                arrows2d!(ax, Point2f(Λ2_inv(h_poisson.times[i]), h_poisson.times[i]), Vec2f(0, -h_poisson.times[i] + arrowspad), minshaftlength=0, color=:orange, alpha=0.5, visible=@lift $visibles2[i])
            end


            axislegend(ax, position=:ct, orientation=:horizontal)
        end

        DOM.div(fig; style="text-align: center")
    end
)

# ╔═╡ 07a7eb40-8b75-45c6-ac8d-79ccddcadaf0
Columns(DOM.div(md"""
!!! warning "Main assumption"
	Particles interact through their empirical distribution
""",
        md"### Main objectives",
        md"""
        1. Find the limit equation [Weiss, 1907]
        """,
        md"""
        2. Prove convergence (rate $n^{-1/2}$) [Sznitman, 1991]
        """,
        md"""
        3. Prove a CLT for the fluctuations [Méléard, 1996]
        """,
        md"""
        4. Derive first-order approximation (rate $n^{-1}$)
        """),
    App() do
        fig = Figure(size=(slideresolution[1] / 3, slideresolution[2]))

        begin
            axleft = Axis(fig[1, 1], title="Particle system", limits=(-3, 3, -3, 3), titlesize=20)
            hidedecorations!(axleft)
            hidespines!(axleft)

            seed!(2)
            positions = rand(MvNormal([0, 0], [1 0; 0 1]), 20)
            scatter!(axleft, positions, markersize=20, color=:orange)
            scatter!(axleft, 0, 0, markersize=20, color=:green)
            for i in 1:size(positions)[2]
                lines!(axleft, [0, positions[1, i]], [0, positions[2, i]], color=:black, alpha=0.1)
            end
        end

        begin
            axright = Axis(fig[2, 1], title="MF approximation", titlesize=20)
            hidedecorations!(axright)
            hidespines!(axright)
            xs = ys = -3:0.05:3
            zs = [pdf(MvNormal([0, 0], [1 0; 0 1]), [x, y]) for x in xs, y in ys]
            heatmap!(axright, xs, ys, zs, colormap=:Oranges_3, alpha=0.7)
            scatter!(axright, [0], [0], markersize=20, color=:green)
            shift = 0
            angles = range(0, 2π, 9)[1:end-1]
            for radius in 0.5:0.5:2.5
                for angle in angles
                    newangle = angle + shift
                    lines!(axright, [0, radius * cos(newangle)], [0, radius * sin(newangle)], color=:black, alpha=0.1)
                end
                shift += angles[2] / 5
            end
        end

        DOM.div(fig)
    end,
    widths=[60, 40])

# ╔═╡ e368d2ee-abe9-40cd-94b1-acaf9ef9482f
DOM.div(
    PlutoUI.LocalResource("assets/raw/shp-raster.png", :height => slideresolution[2] / 3),
    PlutoUI.LocalResource("assets/raw/snfe-raster.png", :height => slideresolution[2] / 3),
    PlutoUI.LocalResource("assets/raw/nfe-raster.png", :height => slideresolution[2] / 3);
    style="text-align: center"
)

# ╔═╡ 795b5f44-3cc6-4f21-96f4-0be0cf015fd9
Columns(DOM.div(@markdown("""
- Neurons are classified into sub-populations
"""),
        @markdown("""
        - Cyclic feedback loop
        - All to all connection
        - Class-wise voltage ``U^k_t``
        - Class-wise intensity ``f_k(U^k_t)``
        """),
        @markdown("""
        - Erlang delay functions:
        ```math
        h_{k\\to k+1}(t) = c_k e^{-\\alpha_k t}\\frac{t^{m_k}}{m_k !}.
        ```
        """
        )),
    DOM.div(PlutoUI.LocalResource("assets/multiclass_v2.png", :width => 0.45 * slideresolution[1]), style="text-align: center"))

# ╔═╡ 66932ace-17fd-4d11-a118-9056c448cd27
App() do
    fig = Figure(size=(slideresolution[1], slideresolution[2] - 350))

    begin # time grid
        kmax = 6
        tmax = 2^3
        ts = tmax * (0:2^kmax) / 2^kmax
    end

    begin # load data
        data = load("assets/coupling_KMT.jld2")
        W = data["Brownian"]
        N = data["Poisson"]
        counting_values = [0.0]
        for i in eachindex(N)
            append!(counting_values, -N[i] .+ [i - 1, i])
        end
        push!(counting_values, -tmax + length(N))
    end

    begin # figure layout parameters
        pad = 0.2
    end

    begin # observables
        btn = Makie.Button(fig[2, 1], label="Start animation", tellwidth=false)

        ts_black = Observable(Float64[])
        ts_gray = Observable(Float64[])
        Ws = Observable(Float64[])
    end

    begin # plots
        ax = Axis(fig[1, 1], limits=((0, tmax) .+ pad .* (-1, 1), (minimum(W), maximum(W)) .+ pad .* (-1, 1)),
            title="KMT coupling", xlabel=L"t", ylabel="",
            xgridvisible=false)
        vlines!(ax, ts_black, color=:black)
        vlines!(ax, ts_gray, color=:gray90)

        lines!(ax, [0; repeat(N, inner=2); tmax], counting_values, alpha=0.7, color=:darkblue)
        initial_W_plot = lines!(ax, ts, W; color=:orange)
        lines!(ax, ts_gray, Ws; color=:orange)
    end

    # Update visualization on each button click
    on(btn.clicks) do clicks

        if clicks == 1
            initial_W_plot.visible = false
            global ks = collect(0:(kmax-1))
            global ids_black = Int[]
            global ids_gray = Int[]
        end

        clicks > 2 * kmax && return
        step = div(clicks, 2)

        if isodd(clicks)  # Odd click: new lines in black
            step == 0 && push!(ids_black, 0, 1)
            for k in ks
                append!(ids_black, 2^k .+ 2^(k - step) * (1:2:2^step))
            end
            popfirst!(ks)

            ts_black[] = ts[ids_black.+1]

            # Update button label
            btn.label[] = step == 0 ? "Compute quantile coupling" : "Compute conditional quantile coupling"
        else  # Even click: turn black lines into gray
            append!(ids_gray, ids_black)
            sort!(ids_gray)
            empty!(ids_black)
            empty!(ts_black[])
            ts_gray[] = ts[ids_gray.+1]
            Ws[] = W[ids_gray.+1]

            # Update button label
            btn.label[] = "Add next ticks"
        end

        clicks == 2 * kmax && (btn.label[] = "Coupling complete!")
    end

    DOM.div(fig; style="text-align: center")
end

# ╔═╡ 74bf35d1-801c-4e01-9705-51578a5a4318
App() do
    fig = Figure(size=(slideresolution[1], slideresolution[2] - 350))

    begin # Load data
        data = load("assets/2CHP-coupling.jld2")

        labels = ["2CHP",
            "diffusion coupled",
            "diffusion independent",
            "MF limit"
        ]

        counting, diffusion_coupled, diffusion_indep, ode = [data[label] for label in labels]
    end

    begin # observables
        gl = GridLayout(fig[1:2, 2], tellheight=false)
        Label(gl[2, 1], "Transparency")
        sl_alpha = Makie.Slider(gl[end+1, 1], range=0:0.01:1, startvalue=1.0, width=100)

        alpha_factor = sl_alpha.value
    end

    begin # axis layout 
        axtop = Axis(fig[1, 1], ylabel=L"U^{1,m'}")
        axbot = Axis(fig[2, 1], ylabel=L"U^{2,m'}", xlabel=L"t")
        linkxaxes!(axtop, axbot)
        hidexdecorations!(axtop, grid=false)
        xlims!(axbot, 0, 140)
    end

    begin # plots
        # 2CHP lines
        for m in 1:4
            lines!(axtop, counting.spike_train[begin:20:end], counting.X[m, begin:20:end];
                color=:darkblue, label="MC-HP", alpha=@lift 1 - $alpha_factor * (m - 1))
        end
        for m in 1:3
            lines!(axbot, counting.spike_train[begin:20:end], counting.X[4+m, begin:20:end];
                color=:darkblue, label="MC-HP", alpha=@lift 1 - $alpha_factor * (m - 1))
        end

        # MF limit lines
        for m in 1:4
            lines!(axtop, ode.time, ode.X[m, :];
                color=:black, label="MF limit", visible=false, alpha=@lift 1 - $alpha_factor * (m - 1))
        end
        for m in 1:3
            lines!(axbot, ode.time, ode.X[4+m, :];
                color=:black, label="MF limit", visible=false, alpha=@lift 1 - $alpha_factor * (m - 1))
        end

        # Independent diffusion lines
        for m in 1:4
            lines!(axtop, diffusion_indep.time, diffusion_indep.X[m, :],
                color=:red, label="Diff. indep.", visible=false, alpha=@lift 1 - $alpha_factor * (m - 1))
        end
        for m in 1:3
            lines!(axbot, diffusion_indep.time, diffusion_indep.X[4+m, :],
                color=:red, label="Diff. indep.", visible=false, alpha=@lift 1 - $alpha_factor * (m - 1))
        end

        # Coupled diffusion lines
        for m in 1:4
            lines!(axtop, diffusion_coupled.time, diffusion_coupled.X[m, :],
                color=:orange, label="Diff. coupled", visible=false, alpha=@lift 1 - $alpha_factor * (m - 1))
        end
        for m in 1:3
            lines!(axbot, diffusion_coupled.time, diffusion_coupled.X[4+m, :],
                color=:orange, label="Diff. coupled", visible=false, alpha=@lift 1 - $alpha_factor * (m - 1))
        end
    end

    # Merge legends from several axes
    # Taken from https://discourse.julialang.org/t/how-to-combine-the-legend-of-multiple-axes-on-same-figure-into-one-in-makie-jl/97299
    plots_in_fig = AbstractPlot[]
    labels_in_fig = AbstractString[]
    for ax in [axtop, axbot]
        pl, lb = Makie.get_labeled_plots(ax, merge=false, unique=false)
        append!(plots_in_fig, pl)
        append!(labels_in_fig, lb)
    end

    ulabels = Base.unique(labels_in_fig)
    mergedplots = [[lp for (i, lp) in enumerate(plots_in_fig) if labels_in_fig[i] == ul]
                   for ul in ulabels]

    Legend(gl[1, 1], mergedplots, ulabels)

    DOM.div(fig; style="text-align: center")
end

# ╔═╡ 7d445019-8982-4293-ae55-9dfa32a7d867
App() do
    fig = Figure(size=(slideresolution[1], slideresolution[2] - 400))

    begin # parameters
        sg = SliderGrid(fig[1, 1],
            (label=L"n", range=50:10:100, startvalue=100),
            (label=L"μ", range=0.0:0.05:1.0, startvalue=0.25),
            (label=L"γ", range=0.0:0.05:1.0, startvalue=0.5),
            (label=L"p", range=0.0:0.05:1.0, startvalue=0.5),
            (label=L"r_+", range=0.0:0.05:1.0, startvalue=0.5),
            tellheight=false)
        n, μ, γ, p, r₊ = [sl.value for sl in sg.sliders]

        warning_text = @lift $μ + $γ > 1 ? "Warning: μ + γ must be less than 1!" : " "
        Label(fig[2, 1], warning_text, tellwidth=false, tellheight=false, color=:red, fontsize=20)
    end

    begin # observables
        model = @lift MarkovChainModel($μ, 1.0 - $γ, $p)
        excitatory = @lift MeanFieldGraph.N2excitatory($n, $r₊)
        simulation = @lift rand($model, $excitatory, 100)

        yticks = @lift [1, floor(Int, $n * $r₊), $n]
    end

    begin # axis and plot
        ax = Axis(fig[:, 2], xlabel=L"t", ylabel="neuron", yticks=yticks, limits=@lift ((nothing, nothing), (1 - 0.05 * $n, 1.05 * $n)))

        poly!((@lift Rect(0.0, 1.0, 100, floor(Int, $n * $r₊) - 0.5)), color=:red, alpha=0.5)
        poly!((@lift Rect(0.0, floor(Int, $n * $r₊) + 0.5, 100, $n - floor(Int, $n * $r₊) - 0.5)), color=:blue, alpha=0.5)
        scatter!(ax, (@lift [p.I for p in findall(transpose($simulation.X))]), markersize=5, color=:black)
    end

    begin # Adjust layout
        colsize!(fig.layout, 1, Fixed(300))
    end

    DOM.div(fig; style="text-align: center")
end

# ╔═╡ 65cb66cf-a7d4-45ef-a2e2-a33f12a2747a
Columns(DOM.div(PlutoUI.LocalResource("assets/colorplot.png", :height => slideresolution[2] / 2),
        PlutoUI.LocalResource("assets/classification_vary_lambda.png", :height => slideresolution[2] / 2);
        style="text-align: center"),
    DOM.div(PlutoUI.LocalResource("assets/classification_vary_N.png", :height => slideresolution[2] / 2),
        PlutoUI.LocalResource("assets/classification_vary_p.png", :height => slideresolution[2] / 2);
        style="text-align: center"))

# ╔═╡ d8e9c413-705a-419b-a3c2-9a2d69d6cfea
DOM.div(PlutoUI.LocalResource("assets/box-outlooks.png", :height => slideresolution[2]), style="text-align: center")

# ╔═╡ dc010f7e-da5d-48e2-a33b-e8305db1cdc8
DOM.div(PlutoUI.LocalResource("assets/droids.jpg", :height => slideresolution[2]);
    style="text-align: center")

# ╔═╡ d64b81f8-855b-4e19-8818-05d10eae80f0
DOM.div(
    PlutoUI.LocalResource("assets/raw/adhp_raster.png", :height => slideresolution[2] / 3),
    PlutoUI.LocalResource("assets/raw/srde_raster.png", :height => slideresolution[2] / 3),
    PlutoUI.LocalResource("assets/raw/rde_raster.png", :height => slideresolution[2] / 3);
    style="text-align: center"
)

# ╔═╡ 0bab1b73-9b42-48a7-9b9f-c6cd738d4992
App() do
    begin # Load data
        data = load("assets/quantiles_abserror.jld2")

        labels = ["q25",
            "q50",
            "q75",
            "q25 large time",
            "q50 large time",
            "q75 large time"
        ]

        q25, q50, q75, q25inf, q50inf, q75inf = [data[label] for label in labels]

        estimatornames = [L"\hat{m}", L"\hat{v}", L"\hat{w}", L"\hat{\mu}", L"\hat{\lambda}", L"\hat{p}"]
    end


    fig = Figure(size=slideresolution)

    # tog = Toggle(fig[2, 2], active=true)
    # visibility = tog.active

    ax = Axis(fig[:, 1], yscale=log10, ylabel="absolute error", xlabel=L"t")
    ylims!(ax, 4.0e-5, 2.7)

    for id in 1:6
        lines!(ax, q50[:, 1], q50[:, id+1]; color=id, colormap=:tab10, colorrange=(1, 10), label=estimatornames[id], visible=id == 1)
        band!(ax, q50[:, 1], q25[:, id+1], q75[:, id+1]; color=id, colormap=:tab10, colorrange=(1, 10), label=estimatornames[id], alpha=0.1, visible=id == 1)
        scatter!(ax, q50[end, 1], q50inf[1, id]; color=id, colormap=:tab10, colorrange=(1, 10), label=estimatornames[id], visible=id == 1)
    end

    Legend(fig[1, 2], ax, merge=true)


    DOM.div(fig; style="text-align: center")
end

# ╔═╡ 78a53f76-521f-4aaf-a798-d186239325f4
# Add color csv style
html"""
<style>
	#gray {
		color: gray;
	}
</style>
<style>
	#red {
		color: red;
	}
</style>
<style>
	#blue {
		color: blue;
	}
</style>
<style>
	#orange {
		color: orange;
	}
</style>
<style>
	#darkblue {
		color: darkblue;
	}
</style>
"""

# ╔═╡ 9a63b59b-1e1d-4b10-8aaa-dc83d7079ba8
# Change fonts
html"""
<style>
	/* Markdown cells */
	pluto-output div {
	    font-size: 25px !important;
	}

	/* Titles */
	pluto-output h1 { font-size: 2.2em !important; font-family: papyrus;}
	pluto-output h2 { font-size: 1.8em !important; }
	pluto-output h3 { font-size: 1.5em !important; }
	pluto-output h4 { font-size: 1.3em !important; }
	pluto-output h5 { font-size: 1.1em !important; }
	pluto-output h6 { font-size: 1em !important; }
</style>
"""

# ╔═╡ 54827aa4-613e-44f4-9fee-1ffb24e72270
RDEplots_text = DOM.div(md"""
$\begin{cases}
	\lambda^{n,i}_t = f(A^{n,i}_{t-}, Ξ^n_t),\\
	A^{n,i}_t := t - \sup\{T \in N^{n,i}, T \leq t\},\\
	Ξ^n_{t} := \frac{1}{n} \sum_{j=1}^{n} \int_{0}^{t-} h(t-s) N^{n,j}(\mathrm{d} s).
\end{cases}
\tag{ADHP}$
""",
    html"</br>",
    html"</br>",
    html"</br>",
    md"""
    $\hspace{-3cm}\begin{cases}
        \{\partial_t + \partial_a\} u^n(t,a) + f(a, \xi^n(t)) u^n(t,a) + \sqrt{\frac{f(a, \xi^n(t)) u^n(t,a)}{n}} W(t,a)= 0,\\
        u^n(t,0) = \int_{0}^{\infty} f(a, \xi^n(t)) u^n(t,a) + \sqrt{\frac{f(a, \xi^n(t)) u^n(t,a)}{n}} W(t,a) \mathrm{d} a,\\
        \xi^n(t) := \int_{0}^{t} h(t-s) u^n(s,0) \mathrm{d} s.
    \end{cases}\tag{SRDE}$
    """,
    html"</br>",
    html"</br>",
    html"</br>",
    md"""
    $\begin{cases}
        \{\partial_t + \partial_a\} u(t,a) + f(a, \xi(t)) u(t,a) = 0,\\
        u(t,0) = \int_{0}^{\infty} f(a, \xi(t)) u(t,a) \mathrm{d} a,\\
        \xi(t) := \int_{0}^{t} h(t-s) u(s,0) \mathrm{d} s.
    \end{cases}
    \tag{RDE}$
    """);

# ╔═╡ 2bdd2fa5-4223-49a5-ace8-1948b8b89f2f
Columns(RDEplots_text,
    DOM.div(PlutoUI.LocalResource("assets/raw/ADHP-hm.png", :height => slideresolution[2] / 3),
        PlutoUI.LocalResource("assets/raw/SRDE-hm.png", :height => slideresolution[2] / 3),
        PlutoUI.LocalResource("assets/raw/RDE-hm.png", :height => slideresolution[2] / 3);
        style="text-align: center"
    )
)

# ╔═╡ 502a2174-437a-4f22-928e-25a6dea73f6b
sensitivity_theorem = md"""
!!! danger "Theorem - C. (2017)"
	If $f(a,\xi) = (\mu + \xi)\mathbf{1}_{[a^*, \infty)}(a)$ then the sensitivity $\partial_\mu \mathrm{act}$ has a unique maximum w.r.t. $\gamma$, denoted by $\gamma_{\rm max}(\mu a^*)$. Furthermore, $\gamma_{\rm max}$ is decreasing and
	
	$\begin{cases}
		\gamma_{\rm max}(\mu a^*) \to 1, & \text{as } \mu a^* \to 0,\\
		\gamma_{\rm max}(\mu a^*) =0, & \text{if } \mu a^* \geq 1/2.
	\end{cases}$
		""";

# ╔═╡ 30ec7e37-701f-43ff-8cdc-2391863b014b
Columns(sensitivity_theorem,
    DOM.div(PlutoUI.LocalResource("assets/RDE-sensitivity_activity.png", :width => 420); style="text-align: center"),
    widths=[60, 40])

# ╔═╡ 3edc08a1-8da5-41ff-bb28-11c1641a2f0f
sensitivity_text = md"""
$\begin{cases}
    \{\partial_t + \partial_a\} u(t,a) + f(a, \xi(t)) u(t,a) = 0,\\
    u(t,0) = \int_{0}^{\infty} f(a, \xi(t)) u(t,a) \mathrm{d} a,\\
    \xi(t) := \int_{0}^{t} h(t-s) u(s,0) \mathrm{d} s.
\end{cases}
\tag{RDE}$

Its stationary version reads as (with $\gamma = \int_0^\infty h(t) \mathrm{d}t$)

$\begin{cases}
    \partial_a u_\infty(a) + f(a, \xi_\infty u_\infty(a)) = 0,\\
    \mathrm{act} := u_\infty(0) = \int_{0}^{\infty} f(a, \xi_\infty u_\infty(a)) \mathrm{d} a,\\
    \xi_\infty := \gamma u_\infty(0).
\end{cases}$

Exponential stability has been studied by Pakdaman et al. (2010), Mischler et al. (2018) and Cáceres et al. (2025).
""";

# ╔═╡ d1328b2d-53a2-4052-8434-65db91077cc6
sensitivity_text

# ╔═╡ 27810bc4-1a6f-49e4-9886-eb6a487857e1
NFEplots_text = DOM.div(md"""
$\begin{cases}
    \text{neuron } i \text{ is at position } x_i \in \mathbb{R}^d,\\ 
    \rho^n(\mathrm{d} x) := n^{-1} \sum_{i=1}^n \delta_{x_i}(\mathrm{d}x),\\
	\lambda^{n,i}_t = f(U^{n,i}_{t-}),\\
	\mathrm{d}U^{n,i}_t = {\color{green}-\alpha U^{n,i}_t} + \frac{1}{n} \sum_{j=1}^{n} w(x_j,x_i) {\color{gold}N^{n,j}(\mathrm{d} t)}.
\end{cases}
\tag{SHP}$
""",
    html"</br>",
    html"</br>",
    md"""
    $\begin{multline}
      \partial_t u^n(t,x) = {\color{green}-\alpha u^n(t,x)} + \int_{\mathbb{R}^d} w(y,x) {\color{gold}f(u^n(t,y))} \rho(\mathrm{d} y) \\
      + \int_{\mathbb{R}^d} w(y,x) {\color{gold}\sqrt{\frac{f(u^n(t,y))}{n}} W(t,y)} \rho(\mathrm{d} y)  
    \end{multline}
    \tag{SNFE}$
    """,
    html"</br>",
    html"</br>",
    html"</br>",
    html"</br>",
    md"""
    $\partial_t u(t,x) = {\color{green}-\alpha u(t,x)} + \int_{\mathbb{R}^d} w(y,x) {\color{gold}f(u(t,y))} \rho(\mathrm{d}y)
    \tag{NFE}$
    """);

# ╔═╡ 745ee1f2-edfd-4982-afef-f910379018bd
Columns(NFEplots_text,
    DOM.div(
        PlutoUI.LocalResource("assets/raw/shp-hm.png", :height => slideresolution[2] / 3),
        PlutoUI.LocalResource("assets/raw/snfe-hm.png", :height => slideresolution[2] / 3),
        PlutoUI.LocalResource("assets/raw/nfe-hm.png", :height => slideresolution[2] / 3);
        style="text-align: center"
    ))

# ╔═╡ 49f682ed-d886-424a-81f2-0220ea44f326
diffusionplots_text = DOM.div(md"""
Neurons are classified into several classes $k$ with cyclic feedback and Erlang delay functions:

$h_{k+1\to k}(t) = c_k e^{-\alpha_k t}\frac{t^{m_k}}{m_k !}.$
""",
    html"</br>",
    md"""
    $\hspace{-0cm}\begin{cases}
        \lambda^{k,i}_t = f_k(U^{k,0}_{t-}),\\
        \mathrm{d} U^{k,m^\prime}_t = \left( -\alpha_k U^{k,m^\prime}_t + U^{k,m^\prime+1}_t \right) \mathrm{d} t, \quad 0\leq m^\prime < m_k,\\
        U^{k,m_k}(\mathrm{d} t) = -\alpha_k U^{k,m_k}_t \mathrm{d} t + c_k \left\{ n_{k+1}^{-1} \sum_{j=1}^{n_{k+1}} N^{k+1, j}(\mathrm{d} t) \right\}.
    \end{cases}
    \tag{MC-HP}$
    """,
    html"</br>",
    html"</br>",
    html"</br>",
    md"""
    $\hspace{-0cm}\begin{cases}
        \mathrm{d} \overline{U}^{k,m^\prime}_t = \left( -\alpha_k \overline{U}^{k,m^\prime}_t + \overline{U}^{k,m^\prime+1}_t \right) \mathrm{d} t, \quad 0\leq m^\prime < m_k,\\
        \mathrm{d}\overline{U}^{k,m_k}_t = -\alpha_k \overline{U}^{k,m_k}_t \mathrm{d} t + c_k\left\{  f_{k+1}(\overline{U}^{k+1,0}_t) \mathrm{d} t + \sqrt{\frac{f_{k+1}(\overline{U}^{k+1,0}_t)}{n_{k+1}}} \mathrm{d} B^{k+1}_t \right\}.
    \end{cases}
    \tag{Diff}$
    """);

# ╔═╡ e9826d8a-6ac5-4f94-b5e7-207b2d6b08f5
diffusionplots_text

# ╔═╡ 432a611c-449f-493c-9d32-8cdff0c1854d
verticalspace(height::Real) = @markdown("""<div style="height: $(string(height))px;"></div>""")

# ╔═╡ ad84ab06-93ea-44d0-9248-4f6b2e204db6
verticalspace(100)

# ╔═╡ b8537acd-e5fc-4be9-87a9-c7eab4458ee4
verticalspace(400)

# ╔═╡ e1e2cc04-faf7-4726-8881-af3eaa5544e6
verticalspace(50)

# ╔═╡ a771033f-f5ec-43cb-bf6e-e74f372002fc
verticalspace(300)

# ╔═╡ a02c4d75-8729-4707-afc4-19d6ae0ae490
verticalspace(200)

# ╔═╡ 810ff3ca-8e10-460f-b7cd-0db6d9c806da
verticalspace(100)

# ╔═╡ f8a4a195-da6e-4b7f-a80a-3c5d5b9d6bd3
verticalspace(200)

# ╔═╡ c8a9699b-3662-4902-99b2-4099e51b7ec4
verticalspace(200)

# ╔═╡ be7baff4-a42a-4bb5-bdea-f12d8b00bbf2
verticalspace(100)

# ╔═╡ ebac01ae-d850-4d93-8499-d4814a84b769
verticalspace(200)

# ╔═╡ 73b77fa8-2e04-471b-97f3-9b0ddace2ffb
verticalspace(100)

# ╔═╡ 2f504d55-91e5-4381-a325-501126841472
verticalspace(100)

# ╔═╡ e02bd698-868c-40b0-9b4f-0e4ceab9b3aa
verticalspace(200)

# ╔═╡ a44b9496-d8b7-4a79-957b-a20f10e92ec3
verticalspace(100)

# ╔═╡ 4907be14-9aa6-42d5-a1c1-ccde39d5ac32
verticalspace(100)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bonito = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
DrWatson = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
GLFW = "f7f18e0c-5ee9-5ccd-a5bf-e8befd85ed98"
MarkdownLiteral = "736d6165-7244-6769-4267-6b50796e6954"
MeanFieldGraph = "f3c20582-38a4-43ca-95fd-b3c95fb7ac78"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PointProcesses = "af0b7596-9bb0-472a-a012-63904f2b4c55"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
Bonito = "~4.1.10"
Distributions = "~0.25.122"
DrWatson = "~2.19.1"
GLFW = "~3.4.5"
MarkdownLiteral = "~0.1.2"
MeanFieldGraph = "~0.2.0"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.73"
PointProcesses = "~0.5.0"
WGLMakie = "~0.13.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "eeb3082c9b4f756706bce8a42dc48c84023fbdd2"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bonito]]
deps = ["Base64", "CodecZlib", "Colors", "Dates", "Deno_jll", "HTTP", "Hyperscript", "LinearAlgebra", "Markdown", "MbedTLS", "MsgPack", "Observables", "OrderedCollections", "Random", "RelocatableFolders", "SHA", "Sockets", "Tables", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "1e8e96a232f4444ee200fc31ef3fc3f8e526411b"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
version = "4.1.10"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChunkCodecCore]]
git-tree-sha1 = "51f4c10ee01bda57371e977931de39ee0f0cdb3e"
uuid = "0b6fb165-00bc-4d37-ab8b-79f91016dbe1"
version = "1.0.0"

[[deps.ChunkCodecLibZlib]]
deps = ["ChunkCodecCore", "Zlib_jll"]
git-tree-sha1 = "cee8104904c53d39eb94fd06cbe60cb5acde7177"
uuid = "4c0bbee4-addc-4d73-81a0-b6caacae83c8"
version = "1.0.0"

[[deps.ChunkCodecLibZstd]]
deps = ["ChunkCodecCore", "Zstd_jll"]
git-tree-sha1 = "34d9873079e4cb3d0c62926a225136824677073f"
uuid = "55437552-ac27-4d47-9aa3-63184e8fd398"
version = "1.0.0"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "e771a63cc8b539eca78c85b0cabd9233d6c8f06f"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "351d6f4eaf273b753001b2de4dffb8279b100769"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.9.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "0037835448781bb46feb39866934e243886d756a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "cb1299fee09da21e65ec88c1ff3a259f8d0b5802"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.4"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "5620ff4ee0084a6ab7097a27ba0c19290200b037"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.4"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cd6756e833c377e0ce9cd63fb97689a255f12323"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "1.33.4+0"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3bc002af51045ca3b47d2e1787d6ce02e68b943a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.122"
weakdeps = ["ChainRulesCore", "DensityInterface", "Test"]

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DrWatson]]
deps = ["Dates", "FileIO", "JLD2", "LibGit2", "MacroTools", "Pkg", "Random", "Requires", "Scratch", "UnPack"]
git-tree-sha1 = "5b6632df14cf24fc2cdb805aab24147001463336"
uuid = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
version = "2.19.1"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "83dc665d0312b41367b7263e8a4d172eac1897f4"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.4"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "eaa040768ea663ca695d442be1bc97edfe6824f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "6.1.3+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "5bfcd42851cf2f1b303f51525a54dc5e98d408a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.15.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW]]
deps = ["GLFW_jll"]
git-tree-sha1 = "40412e58ec374029de3d4ad7c13e1a52aa1e149f"
uuid = "f7f18e0c-5ee9-5ccd-a5bf-e8befd85ed98"
version = "3.4.5"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "1828eb7275491981fa5f1752a5e126e8f26f8741"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.17"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "27299071cc29e409488ada41ec7643e0ab19091f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.17+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "50c11ffab2a3d50192a228c313f05b5b5dc5acb2"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.0+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "815e74f416953c348c9da1d1bc977bbc97c84e18"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.0"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["ChunkCodecLibZlib", "ChunkCodecLibZstd", "FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues"]
git-tree-sha1 = "da2e9b4d1abbebdcca0aa68afa0aa272102baad7"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.6.2"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "368542cde25d381e44d84c3c4209764f05f4ef19"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.6"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MarkdownLiteral]]
deps = ["CommonMark", "HypertextLiteral"]
git-tree-sha1 = "f7d73634acd573bf3489df1ee0d270a5d6d3a7a3"
uuid = "736d6165-7244-6769-4267-6b50796e6954"
version = "0.1.2"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "a370fef694c109e1950836176ed0d5eabbb65479"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.6"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MeanFieldGraph]]
deps = ["Clustering", "Distributions", "LinearAlgebra", "Plots"]
git-tree-sha1 = "d86693fc08717ab94f60f6528f3ac62694e71540"
uuid = "f3c20582-38a4-43ca-95fd-b3c95fb7ac78"
version = "0.2.0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "386b47442468acfb1add94bf2d85365dea10cbab"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f19301ae653233bc88b1810ae908194f07f8db9d"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d922b4d80d1e12c658da7785e754f4796cc1d60d"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.36"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1f7f9bbd5f7a2e5a9f7d96e51c9754454ea7f60b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.4+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "12ce661880f8e309569074a61d3767e5756a199f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.1"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3faff84e6f97a7f18e0dd24373daa229fd358db5"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.73"

[[deps.PointProcesses]]
deps = ["DensityInterface", "Distributions", "LinearAlgebra", "Random", "Statistics", "StatsAPI"]
git-tree-sha1 = "5494ecf679038a212def2d1889e50b5bfeae1b4e"
uuid = "af0b7596-9bb0-472a-a012-63904f2b4c55"
version = "0.5.0"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "34f7e5d2861083ec7596af8b8c092531facf2192"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "8f528b0851b5b7025032818eb5abbeb8a736f853"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "064b532283c97daae49e544bb9cb413c26511f8c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.8"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "8ad2e38cbb812e29348719cc63580ec1dfeb9de4"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.1"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "cec2df8cf14e0844a8c4d770d12347fda5931d72"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.25.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.WGLMakie]]
deps = ["Bonito", "Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "LinearAlgebra", "Makie", "Observables", "PNGFiles", "PrecompileTools", "RelocatableFolders", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "2e3a387f0f71ffb9b1cf5dd48e81581010d667dd"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.13.6"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
"""

# ╔═╡ Cell order:
# ╟─4549be08-558f-4773-b9c4-6fee5a60c97d
# ╟─cdc9af88-ddb8-459e-9041-32d9d02a8e42
# ╟─a10cde06-6230-4d31-8092-0ad0b4e33a0f
# ╠═82e5f7b8-6e9c-4ff0-a49a-22daaf8300cf
# ╠═7a9f4d5f-76a4-4484-b0b4-661ceed69163
# ╟─49e5bea3-ee26-4cd1-b0d1-707b2a5ad337
# ╟─673768b4-e18e-4e68-bd26-315ef9588e08
# ╟─174d611c-1184-426f-8a69-4f9cca0cf989
# ╟─34f58220-87c3-43d0-bca7-a7afbd6cd9c6
# ╟─cbcef609-0781-4728-8a6a-937f82b67c1a
# ╟─84e32427-cc6b-4b35-8a67-a243d8ec8add
# ╠═ad84ab06-93ea-44d0-9248-4f6b2e204db6
# ╠═dd29b1c5-6bd9-48c6-aaa8-de2e951c44ce
# ╟─f04d0ada-3fd5-4441-8d07-67dec5da2504
# ╟─9bdb9101-ec7d-4065-8a17-27e5676951de
# ╠═b8537acd-e5fc-4be9-87a9-c7eab4458ee4
# ╠═8a68c356-c1b7-492b-bdc7-b7ac2dd96a3a
# ╟─0ca0d6a0-2bfd-4f29-aeff-a81b87e3be3b
# ╟─bbbfcb51-81db-4b1d-8382-44c1658658aa
# ╟─e1e2cc04-faf7-4726-8881-af3eaa5544e6
# ╟─78de54a6-6f9a-4a5c-8121-48c9c2a80020
# ╠═a771033f-f5ec-43cb-bf6e-e74f372002fc
# ╠═5a5d3447-fbcc-42a1-92f1-f6bbb6825d3f
# ╟─6007b1d3-17b0-47e8-ab83-6587a66d17e8
# ╟─655cddf9-34f5-41c5-88b8-d095c3dbce5a
# ╟─1dfc7065-05bc-4187-9677-52ecf9a33c1e
# ╟─4e8e58b5-281f-42ab-b3d1-0725ebb5a34f
# ╟─c49ac056-ab65-44ea-b37c-2e9d6426399e
# ╟─68807259-fbb3-4c10-b7b4-5980a59c57e2
# ╠═a02c4d75-8729-4707-afc4-19d6ae0ae490
# ╠═58e967cc-6937-4941-953b-f79c822559c4
# ╟─c448e7c4-34c9-43bf-9171-0f9a1915e8fb
# ╟─dfcf884d-63cd-4e99-a81b-8338760fb5d5
# ╟─357a9f03-912a-4ea0-9d0c-39ff734ebf01
# ╠═810ff3ca-8e10-460f-b7cd-0db6d9c806da
# ╠═7da9719f-2617-4eee-86fa-20dba7728dd2
# ╟─a933cf6b-306f-4940-a4ca-bf493376c515
# ╟─07f4c0df-0d77-4794-bbcd-9a9533c47283
# ╟─661246fe-aea7-48c9-a49d-aed212b10a32
# ╠═f8a4a195-da6e-4b7f-a80a-3c5d5b9d6bd3
# ╠═3f5dba59-d773-45ac-b7e7-f555374a6340
# ╟─9ff8a2e5-4e31-4971-bbd6-fc17d9f1b737
# ╟─0915e18b-0f4b-4fdd-8ab2-ceb26a0eff7b
# ╟─8d01edc8-71ee-4431-a68b-3a3745428176
# ╟─da872269-de92-4003-a2f5-4bb4a79ea228
# ╟─ff5208f3-7116-4dbe-945c-b0ab66e74f9d
# ╟─0cfd87f6-3e2f-423a-894a-7a3223b493c8
# ╟─40dc64c1-cdfc-49a2-aa79-63aea3e52b3c
# ╟─dc8f4bd1-77ea-49c9-8dac-54f8f7064d48
# ╟─07a7eb40-8b75-45c6-ac8d-79ccddcadaf0
# ╠═c8a9699b-3662-4902-99b2-4099e51b7ec4
# ╠═4394c62c-57bf-486a-af2d-cc9c16715e26
# ╟─eeb75f07-3dbe-4b48-8498-e532f2418ed6
# ╟─a3a74379-384e-4e04-9a69-dbaad58a7cd2
# ╟─745ee1f2-edfd-4982-afef-f910379018bd
# ╠═be7baff4-a42a-4bb5-bdea-f12d8b00bbf2
# ╠═ee2d968a-a662-4da4-a17c-19839dffbfc2
# ╟─e368d2ee-abe9-40cd-94b1-acaf9ef9482f
# ╟─05f5212c-44a8-40c7-89b8-ad1abf863fc6
# ╟─e0e76e79-a3cf-4fc8-841b-67a242da0542
# ╟─abe2e182-1f82-46ea-8e09-b845e0671878
# ╟─488ffd0d-cfe6-4cbd-9cd5-7a554cc33868
# ╟─dab4e607-41e7-427f-a9a1-a8b43ca6da8e
# ╟─c4943fad-6567-4b21-a4b2-af9dd2682ba9
# ╟─6d4ed1ad-37e5-4850-b9a3-8e783b867a16
# ╟─cda8bffd-793d-4467-bf02-214456982b91
# ╟─bf5f5043-da32-4136-b536-0ba4ace19703
# ╟─0f9ab402-3e7d-4367-a1e3-ead7968fdb03
# ╟─34950b24-175e-4f66-a220-fcb7df592f36
# ╟─5012cdc0-2399-4ee0-8819-14c75015d9a2
# ╟─795b5f44-3cc6-4f21-96f4-0be0cf015fd9
# ╟─bc594f49-72f5-43a9-a5c3-a8694c0b4b89
# ╠═ebac01ae-d850-4d93-8499-d4814a84b769
# ╠═7d3e5eac-cea0-47f7-840e-485dccc715ec
# ╟─2435405c-8428-4af0-8d4e-4ae099cf3f13
# ╟─fc24f7ef-d9c5-4c52-ae27-266850ce49a2
# ╟─66932ace-17fd-4d11-a118-9056c448cd27
# ╠═73b77fa8-2e04-471b-97f3-9b0ddace2ffb
# ╠═e901002e-e3c7-4d38-be78-c59032d604b6
# ╟─b7542278-c465-4096-921d-aee94dc3adbf
# ╟─0e8bedd2-36e0-4acf-817e-b42a2fd97bbd
# ╟─74bf35d1-801c-4e01-9705-51578a5a4318
# ╟─a4551abd-2aa2-47c3-aea4-d23e81af58d8
# ╟─bc1554fd-6c3c-4217-9df4-1026d5fe9925
# ╟─b38816f9-6cc0-4c16-9d56-cf8f6c2f1ace
# ╟─89b224b0-fe7d-44d5-a641-bcc4a7a472c3
# ╟─533a2768-2c83-483f-9d0d-4dfba4acf8b5
# ╟─5abc34bd-5950-4804-8125-be48d3de4449
# ╟─1b48f36d-fbd2-497d-a175-43aba5f5a6ff
# ╟─51124312-b5e5-4324-92fa-cc00853ea2d6
# ╟─469649cd-a6e5-4d69-8ada-20b566bd3cf1
# ╟─7d445019-8982-4293-ae55-9dfa32a7d867
# ╠═c60f7a22-703e-4cbe-95d8-09039d080d04
# ╟─c562f965-72e2-4a7b-a649-53ac2d83908a
# ╟─31f9380c-2888-4d75-b72b-a1bbbf7e3b06
# ╟─2fc9ed70-0ee5-42ed-b45a-a6acd3e67e91
# ╠═2f504d55-91e5-4381-a325-501126841472
# ╠═2901a2cb-a2c8-48d9-873e-1a6213d29cad
# ╟─00139eee-0abb-4409-8c96-cd44a6395d79
# ╟─32e22997-52f9-46e5-8bdc-2b535953a3c6
# ╟─d3f7f844-14fd-407d-8fd8-94f03afffc2b
# ╠═e02bd698-868c-40b0-9b4f-0e4ceab9b3aa
# ╠═f5f6c19c-b4d1-40b5-8598-fcc66a90dccc
# ╟─8d4f6147-758f-4056-9ce2-3baa1513ee9f
# ╟─65cb66cf-a7d4-45ef-a2e2-a33f12a2747a
# ╠═a44b9496-d8b7-4a79-957b-a20f10e92ec3
# ╠═3188c2b6-f08d-4df2-b2cc-8613cb30145b
# ╟─0d09518c-8e18-45b2-883d-2324c93e68f2
# ╟─42248a46-895a-48e9-8fae-737eb8b63d1d
# ╟─7b101707-4364-4328-af54-3cea58698ab4
# ╟─d8e9c413-705a-419b-a3c2-9a2d69d6cfea
# ╠═4907be14-9aa6-42d5-a1c1-ccde39d5ac32
# ╟─0c1d8ea0-c116-4392-a0ba-7c11a563a036
# ╟─200ae1c8-1521-4f35-98b3-6981d8b7a0d1
# ╟─941709a1-2c58-4d29-b128-84212cee668b
# ╟─dc010f7e-da5d-48e2-a33b-e8305db1cdc8
# ╟─e65af828-655c-4d4e-a040-4599ee179e89
# ╟─eb5975a0-7732-4033-bc66-e7e429611baf
# ╟─5f529f1a-cc46-4662-9fbf-0881151e557f
# ╟─2bdd2fa5-4223-49a5-ace8-1948b8b89f2f
# ╟─a2b2253b-89c1-4bc5-a816-ace56a5099f4
# ╟─d64b81f8-855b-4e19-8818-05d10eae80f0
# ╟─301274ed-109a-44d7-a869-97a49f2bb30b
# ╟─c7b7ff57-5671-447d-9f15-78236f9ffa1a
# ╟─6273cbc1-ca15-43ad-bd66-a7d89232201b
# ╟─d1328b2d-53a2-4052-8434-65db91077cc6
# ╟─30ec7e37-701f-43ff-8cdc-2391863b014b
# ╟─b3275370-4e8c-41fe-bd70-f6bc7eb04a5b
# ╟─be991ec2-e0e5-40f5-8451-c8a1a850e605
# ╟─e9826d8a-6ac5-4f94-b5e7-207b2d6b08f5
# ╟─c019a231-464e-4a2c-95cf-f806cbea5c4a
# ╟─0bab1b73-9b42-48a7-9b9f-c6cd738d4992
# ╟─1eca6154-14ce-4a73-8983-e686997b1b55
# ╟─7ecbde43-ef95-4ed9-942d-149a346ca781
# ╟─5a917ed1-cfc1-40f7-838e-809bdd5a4db3
# ╟─09b7fe1f-4994-42a9-9db0-5f2472262124
# ╟─8062da1a-5b93-4510-918b-dfc8ca767b93
# ╟─28ac8b5f-06c3-4923-b142-24fd7cf8685b
# ╟─932dfc39-958a-4df7-acd0-bcebe912b48e
# ╟─37c5f7e7-5307-4493-8c58-7fce21c969bd
# ╟─0b48d6b6-e4a0-4355-9722-ce4ccd30a560
# ╟─8799c616-5967-49c8-b62a-1aef442bb187
# ╟─41b02bcf-40cf-4234-9bdc-498a29ae6206
# ╟─355cacab-cbed-49a2-8776-9ea4b341c01d
# ╟─ea823dca-9737-4237-911b-d8ba9debf830
# ╟─c71a8f79-762a-4568-93f0-ea1a186eaec8
# ╟─d7a03d41-86af-45fc-895d-b1753d5f2a35
# ╟─dd3c0f63-8bf7-41d2-a958-a681894640ea
# ╟─78a53f76-521f-4aaf-a798-d186239325f4
# ╠═9a63b59b-1e1d-4b10-8aaa-dc83d7079ba8
# ╟─54827aa4-613e-44f4-9fee-1ffb24e72270
# ╟─502a2174-437a-4f22-928e-25a6dea73f6b
# ╟─3edc08a1-8da5-41ff-bb28-11c1641a2f0f
# ╟─27810bc4-1a6f-49e4-9886-eb6a487857e1
# ╟─49f682ed-d886-424a-81f2-0220ea44f326
# ╟─432a611c-449f-493c-9d32-8cdff0c1854d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
