      program main
      
      implicit none
      include "mpp_interface.i"

      character mixture*10, state_model*20, species*10
      integer  i, j, k, ne, ns, nr, var, nm, neq, lin
      parameter(lin=5)
      real(kind=8) T, Tv, P, rho, n, cp, cv, h, mw, e, gamma, a
      real(kind=8) species_x, wdot, species_y, mwi, species_dens, Temps
      dimension species_x(5), wdot(5), species_y(5), mwi(5), Temps(2)
      dimension species_dens(5), energies(10)

c
c     The mixture xml file contains all of the information about the elements,
c     species, reaction mechanism, thermodynamic database and transport models.
c     The state_model controls the chemical and thermodynamic model being used
c     i.e. Equil, ChemNoneq1T, ChemNoneqTTv etc.
      mixture     = "air5"
      state_model = "ChemNonEqTTv"
c
c     Initialise the Mutation++ functions with the mixture file and state model
c     defined above.
      call mpp_initialize(mixture, state_model)
c
c     Can now get information on the number of elements, species, equations etc.
      ne = mpp_nelements()
      ns = mpp_nspecies()
      nr = mpp_nreactions()
      neq = mpp_n_energy_eqns()
      nm = mpp_n_mass_eqns()
c
c     Write the output file header.
      write(*,*) "ne: ",ne,"ns: ",ns,"nr: ",nr,"nen: ",neq,"nm: ",nm      
      write(*,'(A14)',advance='no') "  T(K)        "
      write(*,'(A14)',advance='no') "  Tv(K)        "
      do j = 1,ns
         call mpp_species_name(j, species)
         write(*,'(A12)',advance='no') "Y_"//species
      enddo
      write(*,'(A12)',advance='no') "P(Pa)       "
      write(*,'(A12)',advance='no') "rho(kg/m^3) "
      write(*,'(A12)',advance='no') "n(1/m^3)    "
      write(*,'(A12)',advance='no') "Mw(kg/mol)  "
      write(*,'(A12)',advance='no') "Cp(J/kg-K)  "
      write(*,'(A12)',advance='no') "Cv(J/kg-K)  "
      write(*,'(A12)',advance='no') "Gamma       "
      write(*,'(A12)',advance='no') "a(m/s)      "
      write(*,'(A12)',advance='no') "h(J/kg)     "
      write(*,'(A12)',advance='no') "e(J/kg)     "
      write(*,*)
      
      P = 101325.0
      rho = 1.25
      var = 1
c     Loop over temperature and compute equilibrium properties
      do i = 1,295
         T = dble(i-1)*50.0 + 300.0
         Tv = 300.0
         Temps(1) = T
         Temps(2) = Tv
c
c     Takes inputs of T and P and returns the species mole fractions at
c     equilibrium conditions.
         call mpp_equilibrium_composition(T, P, species_x)
c
c     Converts mole fractions to mass fractions
         call mpp_convert_x_to_y(species_x, species_y)
c         
c     A for loop to give the species densities:
         do k=1,ns
            species_dens(k) = species_y(k)*rho
         enddo
c         
c     Set the state of the mixture using the above species densities and
c     temperature. Once the state of the mixture is set all of the
c     mixture properties can be obtained from the mpp_ functions.
         call mpp_set_state(species_dens, Temps, var)
         
         n   = mpp_number_density()
         rho = mpp_density()
         mw  = mpp_mixture_mw()
         cp  = mpp_mixture_frozen_cp_mass()
         cv  = mpp_mixture_frozen_cv_mass()
         gamma = mpp_mixture_frozen_gamma()
         a = mpp_mixture_frozen_sound_speed()
         h   = mpp_mixture_h_mass()
         e   = mpp_mixture_e_mass()
c
c     Write out properties to output
         write(*,'(2E12.4)',advance='no') T, Tv
         do j = 1,ns
            write(*,'(E12.4)',advance='no') species_y(j)
         enddo
         write(*,'(10E12.4)') P, rho, n, mw, cp, cv, gamma, a, h, e

         call mpp_net_production_rates(wdot)
         do j = 1,ns
            write(*,'(E12.4)') wdot(j)
         enddo
      enddo
c      
c     Clean up the memory stored in the mutation++ library
      call mpp_destroy()

      end program main
