# Derivations

# Semi-Grand Canonical Ensemble

## Monatomic Binary Blend CFC

The system is composed of \f$n_A=n_{A,full} + \lambda\f$ monatomic A molecules
and \f$n_B = n_{B,full} + (1-\lambda)\f$ monatomic B molecules where \f$\lambda
\in (0,1]\f$. 
Microscopic center densities are defined as
\f[
    \begin{align*}
    \hat{\rho}_A(\mathbf r) &= \sum_{j=1}^{n_{A,\textrm{full}}}
        \delta(\mathbf{r} - \mathbf{r}_{j})
        + \lambda \delta(\mathbf{r}-\mathbf{r}_{A,\textrm{part}}) \\
    \hat{\rho}_B(\mathbf r) &= \sum_{j=1}^{n_{B,\textrm{full}}}
        \delta(\mathbf{r} - \mathbf{r}_{j})
        + (1 - \lambda) \delta(\mathbf{r}-\mathbf{r}_{B,\textrm{part}}) \\
    \hat{\rho}_{\textrm{tot}}(\mathbf r) &= \hat{\rho}_A(\mathbf r)
        + \hat{\rho}_B(\mathbf r)
    \end{align*}
\f]

The mean density \f$\rho_0 = n_{tot}/V\f$ remains constant by fixing \f$n_{tot}=n_A +
n_B\f$ and \f$V\f$.
The Flory-like potential energy of interactions between A and B molecules is
given by 
\f[
    \begin{align*}
    \beta U_{1} = 
        \int d \mathbf{r} \int d \mathbf{r}^\prime
        \hat{\rho}_A(\mathbf{r})
        u_{\textrm{eff,1}} (\mathbf{r} - \mathbf{r}^\prime)
        \hat{\rho}_B(\mathbf{r}^\prime)
    \end{align*}
\f]
where
\f[
    \begin{align*}
    u_{\textrm{eff,1}} (\mathbf{r}) = \frac{\chi_{AB}}{\rho_0}
        (h * h)(\mathbf{r})
    \end{align*}
\f]
Local density fluctuations are penalized with a Helfand potential given by
\f[
    \begin{align*}
    \beta U_{2} = 
        \int d \mathbf{r} \int d \mathbf{r}^\prime
        (\hat{\rho}_{\textrm{tot}}(\mathbf{r}) - \rho_0)
        u_{\textrm{eff,2}} (\mathbf{r} - \mathbf{r}^\prime)
        (\hat{\rho}_{\textrm{tot}}(\mathbf{r}^\prime) - \rho_0)
    \end{align*}
\f]
where
\f[
    \begin{align*}
    u_{\textrm{eff,2}} (\mathbf{r}) = \frac{\kappa}{2 \rho_0}
        (h * h)(\mathbf{r})
    \end{align*}
\f]
The semi-grand canonical partition function is given by
\f[
    \begin{align*}
    \Psi = z_1 \int d n_A \int d n_B \int d \lambda
        \frac{\exp[\beta (\mu_A n_A + \mu_B n_B - U_1 - U_2)]}{n_A!n_B!}
        \delta (n_{\textrm{tot}} - n_A - n_B)
    \end{align*}
\f]
We will need a continuous approximation for the factorials in order to
differentiate them with respect to \f$\lambda\f$.
We can use Sterling's approximation which states that
\f[
    \begin{align*}
    n! \approx e^{n \ln n - n}
    \end{align*}
\f]
which lets us rewrite the partition function as
\f[
    \begin{align*}
    \Psi &= z_1 \int d n_A \int d n_B \int d \lambda
        \exp \left[
        \beta (\mu_A n_A + \mu_B n_B - U_1 - U_2)
        - (n_A \ln n_A - n_A) - (n_B \ln n_B - n_B)
        \right]
        \delta (n_{\textrm{tot}} - n_A - n_B)
    \end{align*}
\f]
We'll explicitly evaluate the delta function to eliminate the integral over
\f$n_B\f$.
For the sake of conciseness, we'll continue to use \f$n_B\f$ where
\f$n_B = n_{\textrm{tot}} - n_A\f$.
This gives us
\f[
    \begin{align*}
    \Psi = z_1 \int d n_A \int d \lambda
        \exp ( - \mathcal H)
    \end{align*}
\f]
where
\f[
    \begin{align*}
    \mathcal{H} = 
        \beta (U_1 + U_2)
        - \beta (\mu_A n_A + \mu_B n_B)
        + (n_A \ln n_A - n_A) + (n_B \ln n_B - n_B)
    \end{align*}
\f]
To get the "force" on lambda, we need to evaluate
\f[
    \begin{align*}
    \frac{\partial \lambda}{\partial t} =
        - \zeta \frac{\partial \mathcal H}{\partial \lambda} + \Theta(t)
    \end{align*}
\f]
where \f$\zeta\f$ determines the "diffusivity" for \f$\lambda\f$ and \f$\Theta\f$ is a
Gaussian random variable with variance of \f$2 \zeta \Delta t\f$
\f[
    \begin{align*}
    \frac{\partial \mathcal H}{\partial \lambda} &=
        \frac{\partial \mathcal \beta U_1}{\partial \lambda}
        + \frac{\partial \mathcal \beta U_2}{\partial \lambda}
        - \beta (\mu_A - \mu_B)
        + \ln n_A + 1 - 1
        - \ln n_B - 1 + 1 \\
    &= \frac{\partial \mathcal \beta U_1}{\partial \lambda}
        + \frac{\partial \mathcal \beta U_2}{\partial \lambda}
        - \beta (\mu_A - \mu_B)
        + \ln \frac{n_A}{n_B}
    \end{align*}
\f]
\f[
    \begin{align*}
    \frac{\partial \beta U_1}{\partial \lambda}
        =& \int d \mathbf{r} \int d \mathbf{r}^\prime
        \delta(\mathbf{r} - \mathbf{r}_{A,part})
        u_{\textrm{eff,1}} (\mathbf{r} - \mathbf{r}^\prime)
        \hat{\rho}_B(\mathbf{r}^\prime) \\
        &- \int d \mathbf{r} \int d \mathbf{r}^\prime
        \hat{\rho}_A(\mathbf{r})
        u_{\textrm{eff,1}} (\mathbf{r} - \mathbf{r}^\prime)
        \delta(\mathbf{r}^\prime - \mathbf{r}_{B,part})
    \end{align*}
\f]
\f[
    \begin{align*}
    \frac{\partial \beta U_2}{\partial \lambda}
        =& \int d \mathbf{r} \int d \mathbf{r}^\prime
        ( \delta(\mathbf{r} - \mathbf{r}_{A,part})
        - \delta(\mathbf{r} - \mathbf{r}_{B,part})
        )
        u_{\textrm{eff,2}} (\mathbf{r} - \mathbf{r}^\prime)
        (\hat{\rho}_{tot}(\mathbf{r}^\prime) - \rho_0) \\
        &+ \int d \mathbf{r} \int d \mathbf{r}^\prime
        (\hat{\rho}_{tot}(\mathbf{r}^\prime) - \rho_0)
        u_{\textrm{eff,2}} (\mathbf{r} - \mathbf{r}^\prime)
        ( \delta(\mathbf{r} - \mathbf{r}_{A,part})
        - \delta(\mathbf{r} - \mathbf{r}_{B,part})
        )
    \end{align*}
\f]

## A and B Homopolymer CFC

The system is composed of \f$n_{A,\textrm{tot}}=n_A + \lambda_A\f$ A homopolymers
and \f$n_{B,\textrm{tot}} = n_B + \lambda_B\f$ B homopolymers with \f$N_A\f$ and
\f$N_B\f$ segments each, where both \f$\lambda_A\f$ and \f$\lambda_B \in
(0,1]\f$. \f$n_A\f$ and \f$n_B\f$ represent the number of full A and B
homopolymers, respectively.

Microscopic center densities are defined as
\f[
    \begin{align*}
    \hat{\rho}_A(\mathbf r) &=
        \sum_{j=1}^{n_A}
        \sum_{k=1}^{N_A}
        \delta(\mathbf{r} - \mathbf{r}_{j,k})
        +
        \sum_{k=1}^{N_A}
        \lambda_A \delta(\mathbf{r}-\mathbf{r}_{{n_A}+1,k}) \\
    \hat{\rho}_B(\mathbf r) &=
        \sum_{j=1}^{n_B}
        \sum_{k=1}^{N_B}
        \delta(\mathbf{r} - \mathbf{r}_{j,k})
        +
        \sum_{k=1}^{N_B}
        \lambda_B \delta(\mathbf{r}-\mathbf{r}_{n_B+1,k}) \\
    \hat{\rho}_{\textrm{tot}}(\mathbf r) &= \hat{\rho}_B(\mathbf r)
        + \hat{\rho}_B(\mathbf r)
    \end{align*}
\f]

The mean density \f$\rho_0 = m_{tot}/V\f$ remains constant by fixing
\f$m_{tot}\f$, the total "mass" of the components, and \f$V\f$, the simulation
volume. \f$m_{\textrm{tot}}\f$ is defined as \f$m_{\textrm{tot}}=n_A N_A +
n_B N_B\f$. When transferring mass of amount \f$\Delta m_{AB}\f$ from A to B
homopolymers, the changes to \f$\lambda_A\f$ and \f$\lambda_B\f$ are given by
\f[
    \Delta \lambda_A N_A =\; - \Delta \lambda_B N_B =\; \Delta m_{AB}
\f]

The Flory-like potential energy of interactions between A and B molecules is
given by 
\f[
    \begin{align*}
    \beta U_{1} = 
        \int d \mathbf{r} \int d \mathbf{r}^\prime
        \hat{\rho}_A(\mathbf{r})
        u_{\textrm{eff,1}} (\mathbf{r} - \mathbf{r}^\prime)
        \hat{\rho}_B(\mathbf{r}^\prime)
    \end{align*}
\f]
where
\f[
    \begin{align*}
    u_{\textrm{eff,1}} (\mathbf{r}) = \frac{\chi_{AB}}{\rho_0}
        (h * h)(\mathbf{r})
    \end{align*}
\f]

Local density fluctuations are penalized with a Helfand potential given by
\f[
    \begin{align*}
    \beta U_{2} = 
        \int d \mathbf{r} \int d \mathbf{r}^\prime
        (\hat{\rho}_{\textrm{tot}}(\mathbf{r}) - \rho_0)
        u_{\textrm{eff,2}} (\mathbf{r} - \mathbf{r}^\prime)
        (\hat{\rho}_{\textrm{tot}}(\mathbf{r}^\prime) - \rho_0)
    \end{align*}
\f]
where
\f[
    \begin{align*}
    u_{\textrm{eff,2}} (\mathbf{r}) = \frac{\kappa}{2 \rho_0}
        (h * h)(\mathbf{r})
    \end{align*}
\f]

The semi-grand canonical partition function is given by
\f[
    \begin{align*}
    \Psi =& \int d n_A \int d n_B \int d \lambda_A \int d \lambda_B
        \frac{1}{n_A! n_B! (\lambda_T^3)^{n_A + \lambda_A + n_B + \lambda_B}}
        \\ &\times
        \int d \mathbf{r}^{(n_A + 1) N_A + (n_B + 1) N_B}
        \\ &\times
        \exp[\beta (
            \mu^\prime_A (n_A + \lambda_A) + \mu^\prime_B (n_B + \lambda_B) - U_1 - U_2)
        ]
        \\ &\times
        \delta (
            m_{\textrm{tot}} - (n_A + \lambda_A) N_A - (n_B + \lambda_B) N_B
        )
    \end{align*}
\f]
We can convert the factorials into a continuous form using Sterling's
approximation, which states that
\f[
    \begin{align*}
    n! \approx e^{n \ln n - n}
    \end{align*}
\f]
This lets us rewrite the partition function as
\f[
    \begin{align*}
    \Psi =& \int d n_A \int d n_B \int d \lambda_A \int d \lambda_B
        \int d \mathbf{r}^{(n_A + 1) N_A + (n_B + 1) N_B}
        \\ &\times
        \exp[
            \beta (
                \mu^\prime_A (n_A + \lambda_A) + \mu^\prime_B (n_B + \lambda_B)
                - U_1 - U_2
            )
        ]
        \\ &\times \exp[ - (n_A \ln n_A - n_A) - (n_B \ln n_B - n_B) ]
        \\ &\times \exp[
            - (n_A+\lambda_A) \ln \lambda_T^3
            - (n_B+\lambda_B) \ln \lambda_T^3
        ]
        \\ &\times
        \delta (
            m_{\textrm{tot}} - (n_A + \lambda_A) N_A - (n_B + \lambda_B) N_B
        )
    \end{align*}
\f]
Now let's combine the \f$\lambda_A\f$ and \f$\lambda_A\f$ terms with the
\f$\ln(\lambda_T^3)\f$ terms into effective chemical potentials \f$\mu_A\f$
and \f$\mu_B\f$ according to
\f[
    \begin{align*}
        \beta \mu_A =&\; \beta \mu^\prime_A - \ln(\lambda_T^3) \\
        \beta \mu_B =&\; \beta \mu^\prime_B - \ln(\lambda_T^3)
    \end{align*}
\f]
This gives us
\f[
    \begin{align*}
    \Psi =& \int d n_A \int d n_B \int d \lambda_A \int d \lambda_B
        \int d \mathbf{r}^{(n_A + 1) N_A + (n_B + 1) N_B}
        \\ &\times
        \exp[
            \beta (
                \mu_A (n_A + \lambda_A) + \mu_B (n_B + \lambda_B) - U_1 - U_2
            )
        ]
        \\ &\times
        \exp[ - (n_A \ln n_A - n_A) - (n_B \ln n_B - n_B) ]
        \\ &\times
        \delta (
            m_{\textrm{tot}} - (n_A + \lambda_A) N_A - (n_B + \lambda_B) N_B
        )
    \end{align*}
\f]
Evalating the \f$\delta\f$ function, we can replace \f$\lambda_B\f$ using the
following equation:
\f[
    \begin{align*}
        \lambda_B =& \frac{m_{tot} - (n_A + \lambda_A) N_A}{N_B} - n_B \\
        =& \frac{m_{tot} - m_A}{N_B} - n_B 
    \end{align*}
\f]
The following relationship will be useful:
\f[
    \frac{\partial \lambda_B N_B}{\partial \lambda_A N_A} = -1
\f]
We will keep in mind that \f$\lambda_B\f$ is now defined by this function of
\f$\lambda_A\f$, but for simplicity, we'll leave the term \f$\lambda_B\f$ in
instead of writing out the whole expression. So now we have
\f[
    \begin{align*}
    \Psi =& \int d n_A \int d n_B \int d \lambda_A
        \int d \mathbf{r}^{(n_A + 1) N_A + (n_B + 1) N_B}
        \\ &\times
        \exp[
            \beta (
                \mu_A (n_A + \lambda_A) + \mu_B (n_B + \lambda_B) - U_1 - U_2
            )
        ]
        \\ &\times
        \exp[ - (n_A \ln n_A - n_A) - (n_B \ln n_B - n_B) ]
    \end{align*}
\f]
This can be rewritten as
\f[
    \Psi = \int d n_A \int d n_B \int d \lambda_A
        \int d \mathbf{r}^{(n_A + 1) N_A + (n_B + 1) N_B}
        \exp(-\mathcal{H})
\f]
with an effective Hamiltonian of
\f[
    \begin{align*}
        \mathcal{H} =\;&
            \beta (
                U_1 + U_2 - \mu_A (n_A + \lambda_A) - \mu_B (n_B + \lambda_B)
            )
            \\ &+ n_A (\ln n_A - 1) + n_B (\ln n_B - 1)
    \end{align*}
\f]
To get the "force" on \f$\lambda_A\f$, we need to evaluate
\f[
    \begin{align*}
    \frac{\partial \lambda_A N_A}{\partial t} =
    - \frac{\partial \lambda_B N_B}{\partial t} =
        - \zeta \frac{\partial \mathcal H}{\partial \lambda_A N_A} + \Theta(t)
    \end{align*}
\f]
where \f$\zeta\f$ determines the "diffusivity" for transfer between A and B and
\f$\Theta\f$ is a Gaussian random variable with variance of
\f$2 \zeta \Delta t\f$
\f[
    \begin{align*}
    \frac{\partial \mathcal H}{\partial \lambda_A} =&
        \frac{\partial \mathcal \beta U_1}{\partial \lambda_A}
        + \frac{\partial \mathcal \beta U_2}{\partial \lambda_A}
        - \beta \mu_A
        - \beta \mu_B \frac{\partial \lambda_B}{\partial \lambda_A}
    \\ \frac{\partial \mathcal H}{\partial \lambda_A N_A} =&
        \frac{\partial \mathcal \beta U_1}{\partial \lambda_A N_A}
        + \frac{\partial \mathcal \beta U_2}{\partial \lambda_A N_A}
        - \frac{\beta \mu_A}{N_A}
        + \frac{\beta \mu_B}{N_B}
    \end{align*}
\f]
\f[
    \begin{align*}
    \frac{\partial \beta U_1}{\partial \lambda_A N_A}
        =& \frac{1}{N_A} \sum_{j=1}^{N_A}
        \int d \mathbf{r} \int d \mathbf{r}^\prime
        \delta(\mathbf{r} - \mathbf{r}_{n_A+1,j})
        u_{\textrm{eff,1}} (\mathbf{r} - \mathbf{r}^\prime)
        \hat{\rho}_B(\mathbf{r}^\prime) \\
        &- \frac{1}{N_B} \sum_{j=1}^{N_B}
        \int d \mathbf{r} \int d \mathbf{r}^\prime
        \hat{\rho}_A(\mathbf{r})
        u_{\textrm{eff,1}} (\mathbf{r} - \mathbf{r}^\prime)
        \delta(\mathbf{r}^\prime - \mathbf{r}_{n_B+1,j})
    \end{align*}
\f]
\f[
    \begin{align*}
    \frac{\partial \beta U_2}{\partial \lambda_A N_A}
        =& \int d \mathbf{r} \int d \mathbf{r}^\prime
        \left(
            \frac{1}{N_A} \sum_{j=1}^{N_A}
            \delta(\mathbf{r} - \mathbf{r}_{n_A+1,j})
            - \frac{1}{N_B} \sum_{j=1}^{N_B}
            \delta(\mathbf{r} - \mathbf{r}_{n_B+1,j})
        \right)
        u_{\textrm{eff,2}} (\mathbf{r} - \mathbf{r}^\prime)
        (\hat{\rho}_{tot}(\mathbf{r}^\prime) - \rho_0) \\
        &+ \int d \mathbf{r} \int d \mathbf{r}^\prime
        (\hat{\rho}_{tot}(\mathbf{r}^\prime) - \rho_0)
        u_{\textrm{eff,2}} (\mathbf{r} - \mathbf{r}^\prime)
        \left(
            \frac{1}{N_A} \sum_{j=1}^{N_A}
            \delta(\mathbf{r} - \mathbf{r}_{n_A+1,j})
            - \frac{1}{N_B} \sum_{j=1}^{N_B}
            \delta(\mathbf{r} - \mathbf{r}_{n_B+1,j})
        \right) \\
    \frac{\partial \beta U_2}{\partial \lambda_A N_A}
        =& 2 \int d \mathbf{r} \int d \mathbf{r}^\prime
        \left(
            \frac{1}{N_A} \sum_{j=1}^{N_A}
            \delta(\mathbf{r} - \mathbf{r}_{n_A+1,j})
            - \frac{1}{N_B} \sum_{j=1}^{N_B}
            \delta(\mathbf{r} - \mathbf{r}_{n_B+1,j})
        \right)
        u_{\textrm{eff,2}} (\mathbf{r} - \mathbf{r}^\prime)
        (\hat{\rho}_{tot}(\mathbf{r}^\prime) - \rho_0) \\
    \end{align*}
\f]

## General Semi-Grand Canonical CFC Derivation

Suppose we have a system of \f$ n_s \f$ types of species. The \f$ \chi \f$ and
\f$ \kappa \f$ parts of the nonbonded energy are given by
\f[
    \begin{align*}
      U_{\chi} =& \sum_{I<J} \frac{\chi_{IJ}}{\rho_0}
        \int d \mathbf{r} \int d \mathbf{r^\prime}
        \tilde{\rho}_I(\mathbf{r})
        ( h_I \ast h_J)(\mathbf{r} - \mathbf{r^\prime})
        \tilde{\rho}_J(\mathbf{r^\prime}) \\
      U_{\kappa} =& \frac{\kappa}{2\rho_0}
        \int d \mathbf{r} \left(
            \sum_I (\tilde{\rho}_I \ast h_I)(\mathbf{r}) - \rho_0
        \right)^2
      \\ =& \frac{\kappa}{2 \rho_0}
        \int d \mathbf{r} \int d \mathbf{r^\prime}
        \sum_{I} \tilde{\rho}_I(\mathbf{r})
        ( h_I \ast h_I)(\mathbf{r} - \mathbf{r^\prime})
        \tilde{\rho}_I(\mathbf{r^\prime})
        + \frac{\kappa}{\rho_0}
        \int d \mathbf{r} \int d \mathbf{r^\prime}
        \sum_{I<J} \tilde{\rho}_I(\mathbf{r})
        ( h_I \ast h_J)(\mathbf{r} - \mathbf{r^\prime})
        \tilde{\rho}_J(\mathbf{r^\prime})
        \\ &- 2 \frac{\kappa}{\rho_0} 
        \int d \mathbf{r}
            \rho_0 \sum_I ( \tilde{\rho}_I \ast h_I)(\mathbf{r})
        \\ &+ \frac{\kappa}{2 \rho_0} \int d \mathbf{r} \rho_0^2
    \end{align*}
\f]
Let \f$ m_1 \f$ and \f$ m_2 \f$ represent the mass of 1 whole molecule of
components 1 and 2, respectively. Each component can be composed of any
combination of species. We want to find \f$ \frac{\partial U}{\partial
\lambda_1 m_1} \f$. Let's look at the last two terms in \f$ U_{\kappa} \f$.
The last term is constant, so that one doesn't contribute to the derivative.
The derivative of the second to last term is
\f[
  \begin{align*}
    \frac{\partial}{\partial \lambda_1 m_1}
    \left(
      - 2 \frac{\kappa}{\rho_0} 
      \int d \mathbf{r}
        \rho_0 \sum_I ( \tilde{\rho}_I \ast h_I)(\mathbf{r})
    \right)
    &= - 2 \kappa
    \int d \mathbf{r} \left(
      \frac{1}{m_1} \sum_{j=1}^{N_1} h_{K_{1,j}}(\mathbf{r})
      - \frac{1}{m_2} \sum_{j=1}^{N_2} h_{K_{2,j}}(\mathbf{r})
    \right)
    \\ &= - 2 \kappa
    \int d \mathbf{r} \left( \frac{m_1}{m_1} - \frac{m_2}{m_2} \right)
    \\ &= 0
  \end{align*}
\f]
So now we have
\f[
    \begin{align*}
    \frac{\partial U_\chi}{\partial \lambda_1 m_1} =&
        \frac{1}{m_1}
        \int d\mathbf{r} \int d\mathbf{r^\prime}
        \sum_{J=1}^{n_s}
        \sum_{i=1}^{N_1}
        \frac{\chi_{K_{1,i}J}}{\rho_0}
        \delta(\mathbf{r} - \mathbf{r}_{p1,i})
        (h_{K_{1,i}} \ast h_J)(\mathbf{r} - \mathbf{r}^\prime)
        \tilde{\rho}_J(\mathbf{r}^\prime)
        \\ &- 
        \frac{1}{m_2}
        \int d\mathbf{r} \int d\mathbf{r^\prime}
        \sum_{J=1}^{n_s}
        \sum_{i=1}^{N_2}
        \frac{\chi_{K_{2,i}J}}{\rho_0}
        \delta(\mathbf{r} - \mathbf{r}_{p2,i})
        (h_{K_{2,i}} \ast h_J)(\mathbf{r} - \mathbf{r}^\prime)
        \tilde{\rho}_J(\mathbf{r}^\prime)
    \end{align*}
\f]
and
\f[
    \begin{align*}
    \frac{\partial U_\kappa}{\partial \lambda_1 m_1} =&
        \frac{1}{m_1}
        \int d\mathbf{r} \int d\mathbf{r^\prime}
        \sum_{J=1}^{n_s}
        \sum_{i=1}^{N_1}
        \frac{\kappa}{\rho_0}
        \delta(\mathbf{r} - \mathbf{r}_{p1,i})
        (h_{K_{1,i}} \ast h_J)(\mathbf{r} - \mathbf{r}^\prime)
        \tilde{\rho}_J(\mathbf{r}^\prime)
        \\ &- 
        \frac{1}{m_2}
        \int d\mathbf{r} \int d\mathbf{r^\prime}
        \sum_{J=1}^{n_s}
        \sum_{i=1}^{N_2}
        \frac{\kappa}{\rho_0}
        \delta(\mathbf{r} - \mathbf{r}_{p2,i})
        (h_{K_{2,i}} \ast h_J)(\mathbf{r} - \mathbf{r}^\prime)
        \tilde{\rho}_J(\mathbf{r}^\prime)
    \end{align*}
\f]
Combining them, if we define \f$U_{\chi \kappa} = U_\chi + U_\kappa\f$, we get
\f[
    \begin{align*}
    \frac{\partial U_{\chi \kappa}}{\partial \lambda_1 m_1} =&
        \frac{1}{m_1}
        \int d\mathbf{r} \int d\mathbf{r^\prime}
        \sum_{J=1}^{n_s}
        \sum_{i=1}^{N_1}
        \frac{\chi_{K_{1,i}J} + \kappa}{\rho_0}
        \delta(\mathbf{r} - \mathbf{r}_{p1,i})
        (h_{K_{1,i}} \ast h_J)(\mathbf{r} - \mathbf{r}^\prime)
        \tilde{\rho}_J(\mathbf{r}^\prime)
        \\ &- 
        \frac{1}{m_2}
        \int d\mathbf{r} \int d\mathbf{r^\prime}
        \sum_{J=1}^{n_s}
        \sum_{i=1}^{N_2}
        \frac{\chi_{K_{2,i}J} + \kappa}{\rho_0}
        \delta(\mathbf{r} - \mathbf{r}_{p2,i})
        (h_{K_{2,i}} \ast h_J)(\mathbf{r} - \mathbf{r}^\prime)
        \tilde{\rho}_J(\mathbf{r}^\prime)
    \end{align*}
\f]
