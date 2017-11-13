      program main
      
      implicit none
      include "mpp_interface.i"

      character mixture*10, state_model*20, species*10
      integer  i, j, k, ne, ns, nr, var, nm, neq, lin, its
      integer out_freq
      parameter(lin=5)
      real(kind=8) T, Tv, P, rho, n, cp, cv, h, mw, e, gamma, a
      real(kind=8) species_x, wdot, species_y, mwi, species_dens, Temps
      real(kind=8) wdot_norm, dt, time, energies, rhoe, rhoeve
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
      write(*,*)
      write(*,'(A14)',advance='no') "  Time(s)     "
      write(*,'(A12)',advance='no') "T(K)        "
      write(*,'(A12)',advance='no') "Tv(K)       "
      write(*,'(A12)',advance='no') "P(Pa)       "
      write(*,'(A12)',advance='no') "rho(kg/m^3) "
      do j = 1,ns
         call mpp_species_name(j, species)
         write(*,'(A12)',advance='no') "Y_"//species
      enddo
      do j = 1,ns
         call mpp_species_name(j, species)
         write(*,'(A12)',advance='no') "wdot_"//species
      enddo
      write(*,'(A12)') "wdot(kg/s)  "
      
      P = 101325.0
      rho = 1.25
      var = 1
      T = 4000.0
      Tv = 4000.0
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
c
c     Loop over the set state and wdot functions until we have reached an equilibrium.
      its = 0
      wdot_norm = 1.e2
      time = 0.0
      dt = 8.e-9
      out_freq = 100
 20   if( wdot_norm.gt.1.e1 ) then
c         n   = mpp_number_density()
         rho = mpp_density()
c         mw  = mpp_mixture_mw()
c         cp  = mpp_mixture_frozen_cp_mass()
c         cv  = mpp_mixture_frozen_cv_mass()
c         gamma = mpp_mixture_frozen_gamma()
c         a = mpp_mixture_frozen_sound_speed()
c         h   = mpp_mixture_h_mass()
c     e   = mpp_mixture_e_mass()
         
         call mpp_net_production_rates(wdot)
         wdot_norm = 0.0
         do j = 1,ns
            wdot_norm = wdot_norm + wdot(j)*wdot(j)
            species_dens(j) = species_dens(j) + wdot(j)*dt
         enddo
         wdot_norm = dsqrt(wdot_norm)
c
c     Write out properties to output
         if (MOD(its, out_freq).le.1e-5) then
            write(*,'(5E12.4)',advance='no') time, T, Tv, P, rho 
            do j = 1,ns
               write(*,'(E12.4)',advance='no') species_y(j)
            enddo
            do j = 1,ns
               write(*,'(E12.4)',advance='no') wdot(j)
            enddo
            write(*, '(E12.4)') wdot_norm
         endif
         time = time + dt
         its = its + 1
         if (its.gt.10000) goto 10
         call mpp_species_e_mass( energies )
         rhoe = 0.0
         rhoeve = 0.0
         do j = 1,ns
            rhoe = rhoe + species_dens(j)*energies(j)
            rhoeve = rhoeve + species_dens(j)*energies(j+ns)
         enddo
         var = 0
         temps(1) = rhoe
         temps(2) = rhoeve
         call mpp_set_state(species_dens, temps, var )
         goto 20
      endif   
c      
c     Clean up the memory stored in the mutation++ library
 10   call mpp_destroy()

      end program main
