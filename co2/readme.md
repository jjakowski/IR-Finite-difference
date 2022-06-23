
Here is the the example of IR frequencies  for  CO2 from DFTB+  from Hessian.
The hessian is caluclated using finite difference  DFTB3 and 3ob parameters directly in DFTB+  with this   command;

Driver = SecondDerivatives{
        Atoms= {1:-1 }
       }

We get the following results:

 Vibrational frequencies (all) in cm-1
 1  :         0.0000  (cm-1)
 2  :         0.0000  (cm-1)
 3  :         0.0000  (cm-1)
 4  :         4.2510  (cm-1)
 5  :         4.2511  (cm-1)
 6  :       545.4247  (cm-1)    # bending (double degenerate)
 7  :       545.4247  (cm-1)
 8  :      1363.0969  (cm-1)    # symmetric stretching  (IR inactive, Raman active)
 9  :      2397.1464  (cm-1)    # assymetric stretching (IR active)

Compare also with results from IR-MD


Additional reference:
Note, that NIST database stores a  input & output  Gaussian files  with  vibrational  calculations for CO2 and other molecules at: 
        https://cccbdb.nist.gov/iofiles2.asp

And here is CO2 at B3LYP/6-31G*:
        https://cccbdb.nist.gov/iofiles/co2/m8b1.out.txt
