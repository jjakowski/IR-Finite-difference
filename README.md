# IR-finite-difference
see Peter Atkins, Ronald Friedman, book "*Molecular Quantum Mechanics, Fourth Edition*", Oxford University Press, New York, 2005,  
 page 368 and Example 10.3 on page 361, chapter 10, "*Molecular rotations and vibrations*"


====
 ## IR intensity: 


Mass weighted coordidnates:

                                                        q_i = sqrt(m_i) *x_i

Then hessian in "q"  becomes a *mass weighted hessian*:

                            H_kj = d2 V/ dq_k*dq_j  =   d2 V/dx_k*dx_j  *(m_k*mj)^(-1/2) 
                            EKin = sum_k ( 1/2* m_k* dx_k/dt)   =    1/2*sum_k [ dq_k/dt  ]
                            Epot =  1/2 * sum_k,j [q_k * d2V/dq_k*dq_j * q_j  ]
                                 =  1/2 * sum_k,j [q_k * H_kj* q_j ] 

Diagonalize  H_kj: 

    [cc,ee]  = diag(H_kj)
  which gives normal mode c in mass weighted   coordinates. We convert it back to cartesian normal mode Q:
  
    Q  = cc ./sqrt(m)         

where  

    size(cc,1) = size(m) =  3*Natoms       

and m is  mass vector.
 
We also get reduced mass for k-th mode  from mass vector as : 

                       mass_reduced_k  =    <cc(k) | m | cc(k) > =  sum_i [  cc(k,i)**2 * m_i] 

Now we want to getIR intensity from  dipole  exapnsion::

dipole moment aas a function of  normal mode coordinate Q (displacement from equilibrium):

                       u = u0  + du/dQ * Q  +  1/2* d2u/dQ2 + ....

then

        IR =  | <vib_1| u | vib_0>  | ^2   

        <vib_1| u |  vib_0 > =  <1|u|0> = 
           =  <1|u0|0>  + <1|du/dQ*Q|0 >  +..... 
           =  u0*<1|0>  + du/dQ* <1|Q|0> +  ....
           =  u0* 0     + du/dQ* <1|Q|0> 
           =              du/dQ* <1|Q|0> 

Now the   harmonic oscillator wavefunction as a function of Q and force constant "k": 

       |n> =   |n(Q) >  =   [ sqrt(k/pi) /(2^n *n!) ]^(1/2)  * exp(-k*0.5*Q**2)      * Hn(sqrt(k)*Q ) 
                  Normalization:   N_n = [ sqrt(k/pi) /(2^n *n!_) ]^(1/2)    
                  and Hermites:
                      H0(x)  = 1 
                      H1(x)  = 2 *x     

ground  and first excited vibrational states are:

      |0>   =  [ sqrt(k/pi)]**(1/2)    * exp(-k*0.5*Q^2)        = (k/pi)^(1/4)           * exp(-k*0.5*Q^2)
      |1>   =  [ sqrt(k/pi) / 2]^(1/2) * exp(-k*0.5*Q^2)* (2*Q) =  sqrt(2)* (k/pi)^(1/4) * exp(-k*0.5*Q^2) * Q
    <1|Q|0> =  < Q*exp(-k*0.5*Q^2) |  Q * |  Q*exp(-k*0.5*Q^2)> * sqrt(2*k/pi)
            =    (w*massreduced)^(-1/2)

        IR =    (du/dQ)^2 * 1/(w*massreduced)       #  where : 
         w = sqrt(k/massereduced)


