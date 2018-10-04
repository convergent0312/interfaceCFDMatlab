program two_gammas
    implicit none


    !VARIABLES
    REAL, parameter :: p4 = 2.0 ,p1 = 1.0 , rho4 = 2.0 , rho1 = 1.0
    integer, parameter :: imax = 21
    integer :: i, j
    REAL,parameter :: xmin= 0.0, xmax= 1.0, gamma2 = 1.4, gamma1 = 4.4, p1_infty = 0.0, p2_infty = 6000
    REAL :: dx, temp
    REAL,dimension(1:imax) :: x
    integer :: x0, timestep
    integer, parameter :: maxtimestep = 19
    real,dimension(1:imax) :: u_at_i, rho_at_i, e_at_i, p_at_i, a_at_i, abs_u_plus_a_at_i,dt_at_i, p_infty, E_INTERNAL
    real :: real_dt,dt_smallest
    real :: dt
    real:: epsilon_var
    real :: CFL
    integer :: tmax
    

    real,dimension(1:3,1:imax) :: USTATE_PLUS

    real,dimension(1:3,1:imax) :: dumb


    real,dimension(1:imax) :: rho_plus, u_plus, p_plus, a_plus, e_plus, e_internal_plus, sumplus
    real :: lamda1plus, lamda2plus, lamda3plus
    real :: lamda1plus_plus, lamda2plus_plus, lamda3plus_plus
    real :: lamda1minus_plus, lamda2minus_plus, lamda3minus_plus
    real :: alpha_plus, beta


    real,dimension(1:3,1:imax) :: USTATE_MINUS
    real,dimension(1:imax) :: rho_minus, u_minus, p_minus, a_minus, e_minus, e_internal_minus, summinus
    real :: lamda1minus, lamda2minus, lamda3minus
    real :: lamda1plus_minus, lamda2plus_minus, lamda3plus_minus
    real :: lamda1minus_minus, lamda2minus_minus, lamda3minus_minus
    real :: alpha_minus


    




    !ARRAYS
    REAL, dimension(1:3,1:imax):: USTATE, USTATE_UPDATE !STATE VECTORS
    REAL, dimension(1:imax) :: phi, phi_next, phi_initial, gamma_at_i    


    ! I PLUS ARRAYS
    REAL, dimension(1:3,1:3) :: lamda_plus_i_plus,lamda_minus_i_plus,Ca_i_plus
    REAL,dimension(1:3,1:3) :: Ca_inverse_plus,S_i_plus,S_inverse_plus,xplus,xplusinv,xminus,xminusinv
    REAL,dimension(1:3,1:3) :: Abarplus_plus,Abarminus_plus    
    REAL, dimension(1:3,1:imax) :: FPLUS

    ! I MINUS ARRAYS

    REAL, dimension(1:3,1:3) :: lamda_plus_i_minus,lamda_minus_i_m,Ca_i_minus,Ca_inverse_minus,S_i_minus,S_inverse_minus
    REAL, dimension(1:3,1:3) :: Abarplus_minus,Abarminus_minus    
    REAL, dimension(1:3,1:imax) :: FMINUS



    
    ! INITIAL ARRAYS
    REAL, dimension(1:imax) :: RHO_INT, VELOCITY_INT, PRESSURE_INT, gamma_int


    ! FINAL ARRAYS
    REAL, dimension(1:imax) :: RHO, VELOCITY, PRESSURE, E



    ! STORAGE FOR MATRIX MULTIPLY APLUSPLUS...

    ! i+1/2
    REAL,dimension(1:3,1:3) :: Abarplus_plus_1,Abarplus_plus_2,Abarplus_plus_3
    REAL,dimension(1:3,1:3) :: Abarminus_plus_1,Abarminus_plus_2,Abarminus_plus_3

    ! i-1/2

    REAL,dimension(1:3,1:3) :: Abarplus_minus_1,Abarplus_minus_2,Abarplus_minus_3
    REAL,dimension(1:3,1:3) :: Abarminus_minus_1,Abarminus_minus_2,Abarminus_minus_3


    ! FINITE DIFFERENCE VARIABLES IN DCU

    REAL :: forward_diffX_PHI, backward_diffX_PHI, uplusDCU, uminusDCU


    !calculate dx
    dx = (xmax-xmin)/(imax-1)



    ! epsilon

    epsilon_var = 0.5*dx


    !populating x-array 
    x(1) = xmin
    do i = 1,imax
        x(i+1) = x(i) + dx
    end do

    ! Cut off points, can figure out how to find x0 later

    x0 = 5


    ! populate Phi 

    phi = 0.5*(1+tanh((x-0.5)/(2*epsilon_var)))

    phi_initial = phi
    
    !CALCULATE GAMMA USING 19b
    do  i = 1,imax
        temp = (phi(i)/(gamma1-1)) + ((1-phi(i))/(gamma2-1))
        gamma_at_i(i) = 1 + 1/temp
    end do


    gamma_int = gamma_at_i
    


    !CALCULATE P_INFTY USING 19C

    p_infty = ((gamma_at_i-1)/(gamma_at_i))*(((phi*gamma1*p1_infty)/(gamma1-1))+(((1-phi)*(gamma2*p2_infty))/(gamma2-1)))

    










    ! INITIAL CONDITION

    !! ZERO OUT USTATE
    USTATE = 0.0

    ! RHO 1st
    USTATE(1,1:x0) = 1.362  
    USTATE(1,x0:imax) = 1-0.999*phi(x0:imax)

   

    ! RHO U 2nd
    USTATE(2,1:x0) = USTATE(1,1:x0)*81.25;
    USTATE(2,x0:imax) = USTATE(1,x0:imax)*0.0;

    ! E 
    USTATE(3,1:x0) = (((2.484*10**4)+gamma_at_i(1:x0)*p_infty(1:x0))/(gamma_at_i(1:x0)-1))+0.5*USTATE(1,1:x0)*81.25**2
    USTATE(3,x0:imax) = (((1)+gamma_at_i(x0:imax)*p_infty(x0:imax))/(gamma_at_i(x0:imax)-1))+0.5*USTATE(1,x0:imax)*0.0**2

    ! E INTERNAL

    E_INTERNAL(1:x0) = (USTATE(3,1:x0)/USTATE(1,1:x0))-0.5*81.25**2
    E_INTERNAL(x0:imax) = (USTATE(3,x0:imax)/USTATE(1,x0:imax))-0.5*0.0**2

    ! PRESURE (using Chaklas Eq 19a)

    PRESSURE_INT(1:x0) = ((gamma_at_i(1:x0)-1)*USTATE(1,1:x0)*E_INTERNAL(1:x0)) -(gamma_at_i(1:x0)*p_infty(1:x0))
    PRESSURE_INT(x0:imax) = ((gamma_at_i(x0:imax)-1)*USTATE(1,x0:imax)*E_INTERNAL(x0:imax))-(gamma_at_i(x0:imax)*p_infty(x0:imax))
   
    RHO_INT = USTATE(1,:)
    VELOCITY_INT = USTATE(2,:)/USTATE(1,:)



   !! ZERO OUT ALL ARRAYS


    !! FOR I+1/2
    lamda_plus_i_plus = 0.0
    lamda_minus_i_plus = 0.0
    FPLUS = 0.0
    

     !! FOR I-1/2
    lamda_plus_i_minus = 0.0
    lamda_minus_i_m = 0.0
    FMINUS = 0.0
    Abarplus_minus_1 = 0.0
    Abarplus_minus_2 = 0.0
    Abarplus_minus_3 = 0.0
    Abarplus_minus   = 0.0

    Abarminus_minus_1 = 0.0
    Abarminus_minus_2 = 0.0
    Abarminus_minus_3 = 0.0
    Abarminus_minus   = 0.0

    S_inverse_minus = 0.0
    S_i_minus = 0.0
    Ca_inverse_minus = 0.0
    Ca_i_minus = 0.0


    USTATE_UPDATE = USTATE

    dt_smallest = 0.0
    real_dt = 0.0
    dt = 0.0



    print*, 'ENTER TMAX'
    read*, tmax

    print*, 'ENTER CFL'
    read*, CFL
    
    
    
    !!! STARTING THE MAIN LOOP !!!!!!!


    
    !  MAIN LOOP
    do timestep = 1,tmax

      
    
        ! PART 1: CALCULATE DT  
        do i = 1,imax

            u_at_i = USTATE_UPDATE(2,:)/USTATE_UPDATE(1,:) !rho u by rho 
            rho_at_i = USTATE_UPDATE(1,:)
            e_at_i = USTATE_UPDATE(3,:)
            p_at_i = (e_at_i-0.5*rho_at_i*u_at_i**2)
            a_at_i = sqrt(gamma_at_i(i)*p_at_i/rho_at_i)
            abs_u_plus_a_at_i= abs(u_at_i+a_at_i)
            dt_at_i = dx/abs_u_plus_a_at_i !this is an array of dt at i
            dt_smallest = MINVAL(dt_at_i) !smallest dt in dt array at i
            real_dt = CFL*dt_smallest
        end do 

        dt = real_dt
    
        
    
    
  
    ! PART 2: DEAL WITH U PLUS HALF
        do  i = 2,imax-1

                        
            USTATE_PLUS(:,i)=0.5*(USTATE(:,i)+USTATE(:,i+1))

            rho_plus(:) = USTATE_PLUS(1,:)
           


            u_plus(:) = USTATE_PLUS(2,:)/USTATE_PLUS(1,:)
            e_plus (:) = USTATE_PLUS(3,:)
            e_internal_plus(:) = (e_plus/rho_plus)-(0.5*u_plus**2)

            




            p_plus(i) = (gamma_at_i(i)-1)*rho_plus(i)*E_internal_plus(i)-(gamma_at_i(i)*p_infty(i))

            sumplus(i) = p_plus(i)+p_infty(i)
            sumplus(i) = abs(sumplus(i))            

            
            a_plus(i) = sqrt(real(gamma_at_i(i)*(sumplus(i))/rho_plus(i)))

                                    

            ! ! FORM LAMDAS
            lamda1plus = u_plus(i)
            lamda2plus = u_plus(i)+a_plus(i)
            lamda3plus = u_plus(i)-a_plus(i)
    
            lamda1plus_plus = 0.5*(lamda1plus+abs(lamda1plus))
            lamda2plus_plus = 0.5*(lamda2plus+abs(lamda2plus))
            lamda3plus_plus = 0.5*(lamda3plus+abs(lamda3plus))
            
            lamda1minus_plus = 0.5*(lamda1plus-abs(lamda1plus))
            lamda2minus_plus = 0.5*(lamda2plus-abs(lamda2plus))
            lamda3minus_plus = 0.5*(lamda3plus-abs(lamda3plus))
            

            ! ! LAMDA PLUS FOR I PLUS HALF 
            lamda_plus_i_plus(1,1) = lamda1plus_plus
            lamda_plus_i_plus(2,2) = lamda2plus_plus
            lamda_plus_i_plus(3,3) = lamda3plus_plus

            
            ! ! LAMDA PLUS FOR I MINUS HALF 
            lamda_minus_i_plus(1,1) = lamda1minus_plus
            lamda_minus_i_plus(2,2) = lamda2minus_plus
            lamda_minus_i_plus(3,3) = lamda3minus_plus
            
            ! FORMING Ca, S for I PLUS HALF 
            ! Form Ca plus
            Ca_i_plus(1,1) = 1.0
            Ca_i_plus(1,2) = 0.0
            
            Ca_i_plus(1,3) = -1/(a_plus(i)**2)
            Ca_i_plus(2,1) = 0.0
            Ca_i_plus(2,2) = rho_plus(i)*a_plus(i)
            Ca_i_plus(2,3) = 1.0
            
            Ca_i_plus(3,1) = 0.0
            Ca_i_plus(3,2) = -rho_plus(i)*a_plus(i)
            Ca_i_plus(3,3) = 1.0
            
            
            
            ! Form S plus     
            beta = gamma_at_i(i)-1
            alpha_plus = (u_plus(i)**2)/2
            S_i_plus(1,1) = 1.0
            S_i_plus(1,2) = 0.0
            S_i_plus(1,3) = 0.0
            S_i_plus(2,1) = -u_plus(i)/rho_plus(i)
            S_i_plus(2,2) = 1.0/rho_plus(i)
            S_i_plus(2,3) = 0.0
            S_i_plus(3,1) = alpha_plus*beta
            S_i_plus(3,2) = -u_plus(i)*beta;
            S_i_plus(3,3) = beta;
            
            
            ! Form Ca inverse plus 
            
            Ca_inverse_plus(1,1) = 1.0
            Ca_inverse_plus(1,2) = 1.0/(2*a_plus(i)**2)
            Ca_inverse_plus(1,3) = 1.0/(2.*a_plus(i)**2)
            
            Ca_inverse_plus(2,1) = 0.0
            
            Ca_inverse_plus(2,2) = 1.0/(2*rho_plus(i)*a_plus(i))
            Ca_inverse_plus(2,3) = -1.0/(2*rho_plus(i)*a_plus(i))
            
            Ca_inverse_plus(3,1) = 0.0;
            Ca_inverse_plus(3,2) = 0.5;
            Ca_inverse_plus(3,3) = 0.5;


            
            ! Form S inverse plus 
            
            
            S_inverse_plus(1,1) = 1.0
            S_inverse_plus(1,2) = 0.0
            S_inverse_plus(1,3) = 0.0
            S_inverse_plus(2,1) = u_plus(i)
            S_inverse_plus(2,2) = rho_plus(i)
            S_inverse_plus(2,3) = 0.0
            S_inverse_plus(3,1) = alpha_plus
            S_inverse_plus(3,2) = rho_plus(i)*u_plus(i)
            S_inverse_plus(3,3) = 1.0/beta



            Abarplus_plus_1 = matmul(S_inverse_plus,Ca_inverse_plus)
            Abarplus_plus_2 = matmul(Abarplus_plus_1,lamda_plus_i_plus)
            Abarplus_plus_3 = matmul(Abarplus_plus_2,Ca_i_plus)
            Abarplus_plus   = matmul(Abarplus_plus_3,S_i_plus)



            Abarminus_plus_1 = matmul(S_inverse_plus,Ca_inverse_plus)
            Abarminus_plus_2 = matmul(Abarminus_plus_1,lamda_minus_i_plus)
            Abarminus_plus_3 = matmul(Abarminus_plus_2,Ca_i_plus)
            Abarminus_plus   = matmul(Abarminus_plus_3,S_i_plus)

            FPLUS(:,i) = matmul(Abarplus_plus,USTATE(:,i))+  matmul(Abarminus_plus,USTATE(:,i+1))

            
        end do


        


        ! ! PART 3: DEAL WITH U MINUS HALF
        
        
        

        do  i = 2,imax-1
            USTATE_MINUS(:,i)=0.5*(USTATE(:,i)+USTATE(:,i-1))
            rho_minus(i) = USTATE_MINUS(1,i)            
            u_minus(i) = USTATE_MINUS(2,i)/USTATE_MINUS(1,i)


            e_minus (:) = USTATE_PLUS(3,:)
            e_internal_minus(:) = (e_minus/rho_minus)-(0.5*u_minus**2)

            
            p_minus(i) = (gamma_at_i(i)-1)*rho_minus(i)*e_internal_minus(i)-(gamma_at_i(i)*p_infty(i))

            summinus(i) = p_minus(i)+p_infty(i)
            summinus(i) = abs(summinus(i))            

            
            a_minus(i) = sqrt(real(gamma_at_i(i)*(summinus(i))/rho_minus(i)))
            

            ! FORM LAMBDA PLUS AND MINUS FOR I MINUS HALF 
            lamda1minus = u_minus(i)
            lamda2minus = u_minus(i)+a_minus(i)
            lamda3minus = u_minus(i)-a_minus(i)
            
            lamda1plus_minus = 0.5*(lamda1minus+abs(lamda1minus))
            lamda2plus_minus = 0.5*(lamda2minus+abs(lamda2minus))
            lamda3plus_minus = 0.5*(lamda3minus+abs(lamda3minus))
            
            lamda1minus_minus = 0.5*(lamda1minus-abs(lamda1minus))
            lamda2minus_minus = 0.5*(lamda2minus-abs(lamda2minus))
            lamda3minus_minus = 0.5*(lamda3minus-abs(lamda3minus))
            
            ! LAMDA PLUS FOR I MINUS HALF 
            lamda_plus_i_minus(1,1) = lamda1plus_minus
            lamda_plus_i_minus(2,2) = lamda2plus_minus
            lamda_plus_i_minus(3,3) = lamda3plus_minus
            
            ! LAMDA MINUS FOR I MINUS HALF 
            lamda_minus_i_m(1,1) = lamda1minus_minus
            lamda_minus_i_m(2,2) = lamda2minus_minus
            lamda_minus_i_m(3,3) = lamda3minus_minus
            
            
            
            ! FORMING Ca, S for I MINUS HALF 
            ! Form Ca minus
            Ca_i_minus(1,1) = 1.0
            Ca_i_minus(1,3) = -1/(a_minus(i)**2)
            Ca_i_minus(2,2) = rho_minus(i)*a_minus(i)
            Ca_i_minus(2,3) = 1.0
            Ca_i_minus(3,2) = -rho_minus(i)*a_minus(i)
            Ca_i_minus(3,3) = 1.0
            ! Form S minus     
            beta = gamma_at_i(i)-1
            alpha_minus = (u_minus(i)**2)/2
            S_i_minus(1,1) = 1.0
            S_i_minus(2,1) = -u_minus(i)/rho_minus(i)
            S_i_minus(2,2) = 1.0/rho_minus(i)
            S_i_minus(3,1) = alpha_minus*beta
            S_i_minus(3,2) = -u_minus(i)*beta
            S_i_minus(3,3) = beta
            
            ! Form Ca inverse minus 
            
            Ca_inverse_minus(1,1) = 1.0
            Ca_inverse_minus(1,2) = 1.0/(2*a_minus(i)**2)
            Ca_inverse_minus(1,3) = 1.0/(2*a_minus(i)**2)
            Ca_inverse_minus(2,2) = 1.0/(2*rho_minus(i)*a_minus(i))
            Ca_inverse_minus(2,3) = -1.0/(2*rho_minus(i)*a_minus(i))
            Ca_inverse_minus(3,2) = 0.5
            Ca_inverse_minus(3,3) = 0.5
            
            ! Form S inverse minus 
            
            
            S_inverse_minus(1,1) = 1.0
            S_inverse_minus(2,1) = u_minus(i)
            S_inverse_minus(2,2) = rho_minus(i)
            S_inverse_minus(3,1) = alpha_minus
            S_inverse_minus(3,2) = rho_minus(i)*u_minus(i)
            S_inverse_minus(3,3) = 1.0/beta    
            
            
            

            Abarplus_minus_1 = MATMUL(S_inverse_minus,Ca_inverse_minus)
            Abarplus_minus_2 = MATMUL(Abarplus_minus_1,lamda_plus_i_minus)
            Abarplus_minus_3 = MATMUL(Abarplus_minus_2,Ca_i_minus)
            Abarplus_minus   = MATMUL(Abarplus_minus_3,S_i_minus)


            Abarminus_minus_1 = MATMUL(S_inverse_minus,Ca_inverse_minus)
            Abarminus_minus_2 = MATMUL(Abarminus_minus_1,lamda_minus_i_m)
            Abarminus_minus_3 = MATMUL(Abarminus_minus_2,Ca_i_minus)
            Abarminus_minus   = MATMUL(Abarminus_minus_3,S_i_minus)

            FMINUS(:,i) = matmul(Abarplus_minus,USTATE(:,i-1))+  matmul(Abarminus_minus,USTATE(:,i))          
      

        end do

        !! PART 4: FINITE DIFFERENCE EQUATION
        do i = 2,imax-1
            USTATE_UPDATE(:,i) = USTATE(:,i)-dt/dx*(FPLUS(:,i) - FPLUS(:,i-1) + FMINUS(:,i+1) -FMINUS(:,i))
        end do
            
        

        ! !! PART 5: SETTING BC and PULL

        ! USTATE_UPDATE(:,imax) = USTATE_UPDATE(:,imax-1)
        ! USTATE_UPDATE(:,1) = USTATE_UPDATE(:,2)
        USTATE = USTATE_UPDATE

        RHO = USTATE(1,:)
        VELOCITY = USTATE(2,:)/RHO
        E = USTATE(3,:)

        E_INTERNAL = (E/RHO)-(0.5*VELOCITY**2);

        ! CALCULATE PRESSURE, 19a

        do i = 1,imax
            PRESSURE(i) = ((gamma_at_i(i)-1)*RHO(i)*E_INTERNAL(i))-(gamma_at_i(i))*p_infty(i);
        end do

        
       !! PART 6: CALCULATING PHI

           
        phi_next = phi;
     
        ! Loop and find Uminus, Uplus in the VELOCITY
     
        do i = 2,imax-2
            uplusDCU = MAXVAL(VELOCITY)
            uminusDCU = MINVAL(VELOCITY)
     
            ! the finite difference expressions
            forward_diffX_PHI = phi(i+1)-phi(i)
            backward_diffX_PHI = phi(i)-phi(i-1)
 
         
            !  main equation            
            phi_next(i)= phi(i)-(dt/dx)*(uplusDCU*backward_diffX_PHI+uminusDCU*forward_diffX_PHI) 
                
    
        end do
     
     
        phi = phi_next;
  
 
      !! REACALCULATE GAMMA AT i 


        do i = 1,imax
            temp = (phi(i)/(gamma1-1)) + ((1-phi(i))/(gamma2-1));
            gamma_at_i(i) = 1 + 1/temp;
        end do
            

        !! UPDATE P INFINITY

        p_infty = ((gamma_at_i-1)/(gamma_at_i))*(((phi*gamma1*p1_infty)/(gamma1-1))+(((1-phi)*(gamma2*p2_infty))/(gamma2-1)))
            

          
        
    
    
    
    end do

    


    open(unit = 10, file = 'rdata.dat')
    do i = 1,imax
        write(10,*)x(i),RHO_INT(i), RHO(i)
    end do

    open(unit = 11, file = 'phidata.dat')
    do i = 1,imax
        write(11,*)x(i),phi_initial(i), phi(i)
    end do

    open(unit = 12, file = 'vdata.dat')
    do i = 1,imax
        write(12,*)x(i),VELOCITY_INT(i), VELOCITY(i)
    end do


    open(unit = 13, file = 'prdata.dat')
    do i = 1,imax
        write(13,*)x(i),PRESSURE_INT(i), PRESSURE(i)
    end do

    open(unit = 14, file = 'gammadata.dat')
    do i = 1,imax
        write(14,*)x(i),gamma_int(i), gamma_at_i(i)
    end do

    


    





    
    
   
   

end program two_gammas











        
        ! ! PART 4: FINITE DIFFERENCE EQUATION 

        ! do i = 1,imax

        !     USTATE_UPDATE(:,i) = USTATE(:,i) - (dt/dx)*(FPLUS(:,i)-FMINUS(:,i))

        ! end do

        ! ! PART 5: SETTING BC AND PULLS

        ! USTATE_UPDATE(:,imax) = USTATE_UPDATE(:,imax-1)
        ! USTATE = USTATE_UPDATE
        ! RHO = USTATE(1,:)
        ! VELOCITY = USTATE(2,:)/RHO



        ! ! PART 6: CALCULATE PHI

           
        ! phi_next = phi;       
        

        ! ! Loop and find Uminus, Uplus in the VELOCITY
        
        ! do i = 2,imax-2
        !     uplusDCU = MAXVAL(VELOCITY)
        !     uminusDCU = MINVAL(VELOCITY)
        
        !     ! the finite difference expressions
        !     forward_diffX_PHI = phi(i+1)-phi(i)
        !     backward_diffX_PHI = phi(i)-phi(i-1)
    
            
        !     ! main equation            
        !     phi_next(i)= phi(i)-(dt/dx)*(uplusDCU*backward_diffX_PHI+uminusDCU*forward_diffX_PHI) 
                   
       
        ! end do
        
        
        ! phi = phi_next;
     
        ! ! PART 7: UPDATE GAMMA

        ! do i = 1,imax
        !     temp = (phi(i)/(gamma1-1)) + ((1-phi(i))/(gamma2-1));
        !     gamma_at_i(i) = 1 + 1/temp;
        ! end do