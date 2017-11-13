      program main
      
      implicit none
      include "mpp_interface.i"

      character mixture*10, state_model*20, species*10
      integer  i, j, k, ne, ns, nr, var, nm, neq, lin
      parameter(lin=5)
      real(kind=8) T, Tv, P, rho, n, cp, cv, h, mw, e, gamma, a
      real(kind=8) species_x, wdot, species_y, mwi, species_dens, Temps
      real(kind=8) energies, rhoe, rhoeve
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
      write(*,'(A12)',advance='no') "Tv(K)        "
      write(*,'(A12)',advance='no') "rho(kg/m^3) "
c      do j = 1,ns
c         call mpp_species_name(j, species)
c         write(*,'(A12)',advance='no') "Y_"//species
c      enddo
      do j = 1,ns
         call mpp_species_name(j, species)
         write(*,'(A12)',advance='no') "e_tr_"//species
      enddo
      do j = 1,ns
         call mpp_species_name(j, species)
         write(*,'(A12)',advance='no') "e_ve_"//species
      enddo
      write(*,'(A12)',advance='no') "rhoe(J/kg)     "
      write(*,'(A12)',advance='no') "rhoeve(J/kg)     "
      write(*,*)
      
      P = 101325.0
      rho = 1.25
      var = 1
c     Loop over temperature and compute equilibrium properties
      T = 2000.0
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
      do j = 1,ns*2
         energies(j) = 0.0d0
      enddo
      call mpp_species_e_mass( energies )
      rhoe = 0.0
      rhoeve = 0.0
      do j = 1,ns
         rhoe = rhoe + species_dens(j)*energies(j)
         rhoeve = rhoeve + species_dens(j)*energies(j+ns)
      enddo
      rho = mpp_density()
c
c     Write out properties to output
      write(*,'(3E12.4)',advance='no') T, Tv, rho
      do j = 1,ns
         write(*,'(E12.4)',advance='no') energies(j) - energies(j+ns)
      enddo
      do j = 1,ns
         write(*,'(E12.4)',advance='no') energies(j+ns)
      enddo
      write(*,'(3E12.4)') rhoe, rhoeve
      write(*,*)
c
c     Set the state using the energies and we should recover the same result.
      var = 0
      temps(1) = rhoe
      temps(2) = rhoeve
      call mpp_set_state(species_dens, temps, var )
      
      call mpp_species_e_mass( energies )
      rhoe = 0.0
      rhoeve = 0.0
      do j = 1,ns
         rhoe = rhoe + species_dens(j)*energies(j)
         rhoeve = rhoeve + species_dens(j)*energies(j+ns)
      enddo
      rho = mpp_density()
c
c     Write out properties to output
      write(*,'(3E12.4)',advance='no') T, Tv, rho
      do j = 1,ns
         write(*,'(E12.4)',advance='no') energies(j) - energies(j+ns)
      enddo
      do j = 1,ns
         write(*,'(E12.4)',advance='no') energies(j+ns)
      enddo
      write(*,'(3E12.4)') rhoe, rhoeve
      write(*,*)
c      
c     Clean up the memory stored in the mutation++ library
      call mpp_destroy()

      end program main
