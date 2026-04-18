#' Simulated Eisenia fetida Growth Data
#'
#' Simulated growth dataset for 21 individuals of the earthworm
#' *Eisenia fetida* measured weekly over 12 weeks under controlled
#' laboratory conditions (20 °C, ad libitum feeding), following a
#' standard OECD 222 test protocol design.  DEB parameters are based on
#' the *E. fetida* entry in the Add-my-Pet collection (AmP;
#' Marques et al., 2018).  Individual variation in \eqn{\{p_{Am}\}} was
#' drawn from a lognormal distribution with CV \eqn{\approx} 10\%, and
#' Gaussian measurement error \eqn{\sigma_L = 0.015} cm was added to
#' the structural length.
#'
#' @format A data frame with 273 rows and 3 columns:
#' \describe{
#'   \item{id}{Individual identifier (integer, 1--21)}
#'   \item{time}{Observation time in days (0, 7, 14, \ldots, 84)}
#'   \item{length}{Structural length in cm (= \eqn{V^{1/3}}), with
#'     measurement error}
#' }
#' @source Simulated from the standard DEB model (Kooijman, 2010,
#'   Eqs. 2.4--2.6) with parameters:
#'   \eqn{\{p_{Am}\} = 5.0} J/d/cm\eqn{^2},
#'   \eqn{[p_M] = 0.5} J/d/cm\eqn{^3},
#'   \eqn{\kappa = 0.75},
#'   \eqn{v = 0.2} cm/d,
#'   \eqn{[E_G] = 400} J/cm\eqn{^3}.
#'   Individual variation: \eqn{\{p_{Am}\}_j \sim
#'   \mathrm{LogNormal}(\log 5.0,\, 0.1)},
#'   \eqn{L_{0,j} \sim \mathrm{LogNormal}(\log 0.1,\, 0.05)}.
#'   Observation error: \eqn{\sigma_L = 0.015} cm.
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press.
#' \doi{10.1017/CBO9780511805400}
#'
#' Marques, G.M., Augustine, S., Lika, K., Pecquerie, L., Domingos, T.
#' and Kooijman, S.A.L.M. (2018). The AmP project: comparing species on
#' the basis of dynamic energy budget parameters. *PLOS Computational
#' Biology*, 14(5), e1006100. \doi{10.1371/journal.pcbi.1006100}
#'
#' OECD (2016). Test No. 222: Earthworm Reproduction Test (*Eisenia
#' fetida* / *Eisenia andrei*). OECD Guidelines for the Testing of
#' Chemicals.
#'
#' @examples
#' data(eisenia_growth)
#' head(eisenia_growth)
"eisenia_growth"

#' Simulated Folsomia candida Reproduction Data
#'
#' Simulated 28-day reproduction dataset for the springtail *Folsomia
#' candida* exposed to 5 cadmium concentrations (0, 10, 50, 100,
#' 500 mg Cd/kg dry soil) with 6 replicates per concentration.  The
#' experimental design follows ISO 11267 (springtail reproduction test).
#' Toxicant effects were simulated using a simple NEC-based model:
#' expected offspring \eqn{= 120 \times \max(1 - b_w \max(C - z_w, 0),
#' \, 0)}, with Poisson sampling around the expected value.
#'
#' @format A data frame with 30 rows and 5 columns:
#' \describe{
#'   \item{id}{Replicate identifier (integer, 1--30)}
#'   \item{concentration}{Cadmium concentration in mg/kg dry soil}
#'   \item{t_start}{Start of observation interval (day 0)}
#'   \item{t_end}{End of observation interval (day 28)}
#'   \item{count}{Number of juveniles produced in the interval}
#' }
#' @source Simulated with NEC \eqn{z_w = 30} mg/kg, effect intensity
#'   \eqn{b_w = 0.005} (mg/kg)\eqn{^{-1}}, and control offspring
#'   count \eqn{\approx 120}.
#'
#' @references
#' ISO 11267 (2014). *Soil quality — Inhibition of reproduction of
#' Collembola (Folsomia candida) by soil contaminants*.
#'
#' @examples
#' data(folsomia_repro)
#' head(folsomia_repro)
"folsomia_repro"

#' Simulated DEBtox Growth Data
#'
#' Simulated growth data under toxicant exposure for 4 concentration
#' groups (0, 20, 80, 200 arbitrary units) with 10 individuals per group
#' measured weekly over 6 weeks.  Designed for fitting the DEBtox (TKTD)
#' model in [bdeb_tox()].  The toxicant acts through stress on
#' assimilation: the effective assimilation rate is reduced by a factor
#' \eqn{\max(1 - b_w \max(C_w - z_w, 0), \, 0)} following the DEBtox
#' framework of Jager et al. (2006).
#'
#' @format A data frame with 280 rows and 4 columns:
#' \describe{
#'   \item{id}{Individual identifier (integer, 1--40)}
#'   \item{concentration}{External toxicant concentration (arbitrary units)}
#'   \item{time}{Observation time in days (0, 7, 14, \ldots, 42)}
#'   \item{length}{Structural length in cm (with measurement error)}
#' }
#' @source Simulated with NEC \eqn{z_w = 15}, effect intensity
#'   \eqn{b_w = 0.003}, and base DEB parameters identical to
#'   [eisenia_growth].
#'
#' @references
#' Jager, T., Heugens, E.H.W. and Kooijman, S.A.L.M. (2006). Making
#' sense of ecotoxicological test results: towards application of
#' process-based models. *Ecotoxicology*, 15(3), 305--314.
#' \doi{10.1007/s10646-006-0060-x}
#'
#' @examples
#' data(debtox_growth)
#' head(debtox_growth)
"debtox_growth"

#' Eisenia fetida Growth Data (Neuhauser 1980)
#'
#' Real experimental growth data for *Eisenia fetida* on activated
#' sludge at 25 °C from Neuhauser, Hartenstein & Kaplan (1980),
#' obtained via the EGrowth database (Mathieu, 2018).  The dataset
#' comprises 37 group-mean body mass measurements (20 worms per
#' replicate) over 250 days, from hatchling (~8 mg) to adult
#' (~2350 mg).  Wet mass was converted to structural length via
#' \eqn{L = (W / d_V)^{1/3}} with tissue density
#' \eqn{d_V = 1.05} g cm\eqn{^{-3}}.
#'
#' @format A data frame with 37 rows and 3 columns:
#' \describe{
#'   \item{id}{Individual/group identifier (always 1)}
#'   \item{time}{Time in days since start (1--250)}
#'   \item{length}{Structural length in cm (0.20--1.31)}
#' }
#'
#' @source EGrowth database, curve ID gr0226.
#'   \url{https://github.com/JeromeMathieuEcology/EGrowth}
#'
#' @references
#' Neuhauser, E.F., Hartenstein, R. and Kaplan, D.L. (1980). Growth
#' of the earthworm *Eisenia foetida* in relation to population density
#' and food rationing. *Oikos*, 35, 93--98.
#' \doi{10.2307/3544730}
#'
#' Mathieu, J. (2018). EGrowth: a global database on intraspecific
#' body growth variability in earthworm. *Soil Biology and
#' Biochemistry*, 122, 71--80. \doi{10.1016/j.soilbio.2018.04.004}
#'
#' @examples
#' data(eisenia_neuhauser)
#' plot(eisenia_neuhauser$time, eisenia_neuhauser$length,
#'      xlab = "Days", ylab = "Structural length (cm)")
"eisenia_neuhauser"

#' Eisenia andrei Cadmium Toxicity Data (Van Gestel 1991)
#'
#' Real experimental growth data for *Eisenia andrei* exposed to
#' cadmium in natural soil at 23 °C from Van Gestel et al. (1991),
#' obtained via the EGrowth database (Mathieu, 2018).  Five
#' concentration groups (0, 10, 32, 100, 320 mg Cd/kg) with 7 time
#' points each over 85 days.  Group-mean body mass was converted to
#' structural length via \eqn{L = (W / d_V)^{1/3}} with
#' \eqn{d_V = 1.05} g cm\eqn{^{-3}}.
#'
#' @format A data frame with 35 rows and 4 columns:
#' \describe{
#'   \item{id}{Concentration group identifier (1--5)}
#'   \item{time}{Time in days (0--85)}
#'   \item{length}{Structural length in cm}
#'   \item{concentration}{Cadmium concentration (mg Cd/kg soil)}
#' }
#'
#' @source EGrowth database, curve IDs gr0119--gr0123.
#'   \url{https://github.com/JeromeMathieuEcology/EGrowth}
#'
#' @references
#' Van Gestel, C.A.M., Van Dis, W.A., Dirven-Van Breemen, E.M.,
#' Sparenburg, P.M. and Baerselman, R. (1991). Influence of cadmium,
#' copper, and pentachlorophenol on growth and sexual development of
#' *Eisenia andrei* (Oligochaeta; Annelida). *Biology and Fertility
#' of Soils*, 12, 117--121. \doi{10.1007/BF00341486}
#'
#' @examples
#' data(eisenia_cd)
#' plot(eisenia_cd$time, eisenia_cd$length,
#'      col = as.factor(eisenia_cd$concentration),
#'      xlab = "Days", ylab = "Structural length (cm)")
"eisenia_cd"
