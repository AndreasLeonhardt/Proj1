#include "trialfct_molecule.h"

TrialFct_molecule::TrialFct_molecule()
{
}
TrialFct_molecule::TrialFct_molecule(Config * parameters) : function(parameters)
{
}






// calculate the value of the trial function the position given in the variable position.
// Slater determinant of hydrogen wave functions and Jastrow factor.
double TrialFct_molecule::getValue(positions * R)
{
        double result = 0.0;

        for(int i=0;i<nParticles;i++)
        {
            result += hydrogen(i,0,R);
        }



// SLATER DETERMINANT FOR SINGLE ATOM
//        // create Slater matrix
//        mat Slaterup  = zeros(nParticleshalf,nParticleshalf);
//        mat Slaterdown  = zeros(nParticleshalf,nParticleshalf);



//        for (int i=0;i<nParticleshalf;i++)
//        {
//            for (int j=0;j<nParticleshalf;j++)
//            {
//                Slaterup(i,j)   = hydrogen(i,j,R);
//                Slaterdown(i,j) = hydrogen(i+nParticleshalf,j,R);
//            }
//        }
//        result = det(Slaterup)*det(Slaterdown);


       //  jastrow factor

        double jastrow = 0.0;
        double dist;

        for (int i=0;i<nParticleshalf;i++)
        {   // equal spin, factor 1/4
            for (int j=0;j<i;j++)
            {
                dist = R->get_rr(i-1,j);
                jastrow += dist/(4+4*funcParameters[1]*dist);

                dist = R->get_rr(nParticleshalf+i-1,nParticleshalf+j);
                jastrow += dist/(4+4*funcParameters[1]*dist);
            }

            // opposite spin, factor 1/2
            for (int j=nParticleshalf;j<nParticles;j++)
            {
                dist = R->get_rr(j-1,i);
                jastrow += dist/(2+2*funcParameters[1]*dist);
            }

        }

        result *= exp(jastrow);

        return result;
}


// not used for hydrogen molecule
void TrialFct_molecule::setSlaterinv(positions * R)
{

//    // create Slater matrix
//    mat Slaterup  = zeros(nParticleshalf,nParticleshalf);
//    mat Slaterdown  = zeros(nParticleshalf,nParticleshalf);



//    for (int i=0;i<nParticleshalf;i++)
//    {
//        for (int j=0;j<nParticleshalf;j++)
//        {
//            Slaterup(i,j)   = hydrogen(i,               j,R);
//            Slaterdown(i,j) = hydrogen(i+nParticleshalf,j,R);
//        }
//    }

//    // calculate the inverse of the Slatermatrix for spin up and down
//    // respectively. This is used for later updates of the position.
//    inverseSlaterDown = inv(Slaterdown);
//    inverseSlaterUp   = inv(Slaterup);

}


// after moving one particle, this updates the inverse of the Slater matrix,
// however, there are some issues, since the updatet matrix differs from the one
// set up from the new position with teh above function setSlaterinv.
// not used for hydrogen molecule
void TrialFct_molecule::updateSlaterinv(int particleNumber, positions* Rnew, double ratio)
{
//    mat newslatermat = zeros(nParticleshalf,nParticleshalf);
//    double sum;
//    // if particle has spin up, just update SlaterUp
//    if(particleNumber<nParticleshalf)
//    {
//        for (int k=0;k<nParticleshalf;k++)
//        {
//            for (int j=0;j<nParticleshalf;j++)
//            {
//                if(particleNumber==j)
//                {
//                    newslatermat(k,j)=inverseSlaterUp(k,j)/ratio;
//                }
//                else
//                {   sum=0.0;
//                    for (int l=0;l<nParticleshalf;l++)
//                    {
//                        sum+= hydrogen(particleNumber,l,Rnew)*inverseSlaterUp(l,j);
//                    }
//                    newslatermat(k,j) =inverseSlaterUp(k,j) -inverseSlaterUp(k,particleNumber)*sum/ratio;
//                }
//            }
//        }
//        inverseSlaterUp = newslatermat;
//    }
//    // if it has spin down, just update spin down
//    else
//    {
//        for (int k=0;k<nParticleshalf;k++)
//        {
//            for (int j=0;j<nParticleshalf;j++)
//            {
//                if(particleNumber==j+nParticleshalf)
//                {
//                    newslatermat(k,j)=inverseSlaterDown(k,j)/ratio;
//                }
//                else
//                {   sum=0.0;
//                    for (int l=0;l<nParticleshalf;l++)
//                    {
//                        sum+= hydrogen(particleNumber,l,Rnew)*inverseSlaterDown(l,j);
//                    }
//                    newslatermat(k,j)=inverseSlaterDown(k,j)-inverseSlaterDown(k,particleNumber-nParticleshalf)/ratio*sum;
//                }
//            }
//        }
//        inverseSlaterDown = newslatermat;
//    }
}


//================================================================================================
// (laplace \Psi)/\Psi
// This includes several parts:
//   (laplace |D|) / |D|
// + div((grad Jastrow)/Jastrow)
// + ((grad Jastrow)/Jastrow)^2
// + 2* (grad |D|)/|D| * (grad Jastrow)/Jastrow

double TrialFct_molecule::getDivGradOverFct(int particleNumber, positions *R)
{
    double result = divgradhydrogen(particleNumber,0,R)/hydrogen(particleNumber,0,R);

    //-------------------------------------------------------------------------------
    // div( grad(Jastrow)/ Jastrow)

    double b = funcParameters[1];
    double rr;
    if(particleNumber<nParticleshalf)
    {
        int i;
        for (i=0;i<particleNumber;i++)
        {
            rr=R->get_rr(particleNumber-1,i);
            result += .5 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }

        for (i=particleNumber+1;i<nParticleshalf;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += .5 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }

        for (i=nParticleshalf;i<nParticles;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += 1 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }
    }

    else
    {
        int i;
        for (int i=0;i<nParticleshalf;i++)
        {
            rr=R->get_rr(particleNumber-1,i);
            result += 1 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }
        for (i=nParticleshalf;i<particleNumber;i++)
        {
            rr=R->get_rr(particleNumber-1,i);
            result += .5 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }

        for (i=particleNumber+1;i<nParticles;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += .5 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }
    }
    //------------------------------------------------------------------

    // (div Jastrow)^2
    vec jastrowgradient = GradJastrow(particleNumber,R);
    result += dot(jastrowgradient,jastrowgradient);

    // 2*(div Jastrow) * (div SlaterDet)
    result += 2*dot(GradSlater(particleNumber,R),jastrowgradient);

    return result;
}








//===========================================================================================
// returns the quantum force. using closed form expressions
// and the inverse slater determinant matrix. Be careful, that this matrix is
// set to the right value.
vec TrialFct_molecule::quantumForce(int particleNumber, positions *R)
{
    vec result = zeros(ndim);

    result = GradSlater(particleNumber,R);
    //cout <<result<<endl;

    result += GradJastrow(particleNumber,R);
   // cout <<result << endl;

    return 2*result;
}

//=====================================================================================================
// calculate the ratio of two Slaterdeterminants, that differ only in one particle
double TrialFct_molecule::SlaterRatio(int particleNumber ,positions * Rold,positions * Rnew)
{
     double result = hydrogen(particleNumber, 0, Rnew);
     result /= hydrogen(particleNumber,0,Rold);
    return result;
}


double TrialFct_molecule::JastrowRatio(int particleNumber, positions * Rold, positions * Rnew)
{
    // jastrow factor new
    double jastrow = 0.0;
    double dist;



    if (particleNumber<nParticleshalf)
    {
        for (int i=0;i<particleNumber;i++)
        {
            dist = Rnew->get_rr(particleNumber-1,i);
            jastrow += dist/(4+4*funcParameters[1]*dist);

            dist = Rold->get_rr(particleNumber-1,i);
            jastrow -= dist/(4+4*funcParameters[1]*dist);
        }
        for (int i=particleNumber+1;i<nParticleshalf;i++)
        {
            dist = Rnew->get_rr(i-1,particleNumber);
            jastrow += dist/(4+4*funcParameters[1]*dist);

            dist = Rold->get_rr(i-1,particleNumber);
            jastrow -= dist/(4+4*funcParameters[1]*dist);
        }
        for (int i=nParticleshalf;i<nParticles;i++)
        {
            dist = Rnew->get_rr(i-1,particleNumber);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = Rold->get_rr(i-1,particleNumber);
            jastrow -= dist/(2+2*funcParameters[1]*dist);
        }
    }
    else
    {
        for (int i=0;i<nParticleshalf;i++)
        {
            dist = Rnew->get_rr(particleNumber-1,i);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = Rold->get_rr(particleNumber-1,i);
            jastrow -= dist/(2+2*funcParameters[1]*dist);
        }
        for (int i=nParticleshalf;i<particleNumber;i++)
        {
            dist = Rnew->get_rr(particleNumber-1,i);
            jastrow += dist/(4+4*funcParameters[1]*dist);

            dist = Rold->get_rr(particleNumber-1,i);
            jastrow -= dist/(4+4*funcParameters[1]*dist);
        }
        for (int i=particleNumber+1;i<nParticles;i++)
        {
            dist = Rnew->get_rr(i-1,particleNumber);
            jastrow += dist/(4+4*funcParameters[1]*dist);

            dist = Rold->get_rr(i-1,particleNumber);
            jastrow -= dist/(4+4*funcParameters[1]*dist);
        }
    }



    return exp(jastrow);
}

// calculates the divergence of the Jastrow factor divided by the jastrow factor
// That is sum_{i!=k} \frac{a}{ (1+\beta*r_{ki})^2 }  \frac{ \vec{r}_{ki} }{ r_{ki} }
// where a is 0.5 for spin_i = spin_k and
//            0.25 for spin_i!=spin_k
vec TrialFct_molecule::GradJastrow(int particleNumber, positions * R)
{
    vec result = zeros(ndim);
    double b = funcParameters[1];
    double rr;
    if(particleNumber<nParticleshalf)
    {
        for (int i=0;i<particleNumber;i++)
            {
            rr=R->get_rr(particleNumber-1,i);
            result += 0.25 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                    /((1+b*rr)*(1+b*rr)*rr);
            }

        for (int i=particleNumber+1;i<nParticleshalf;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += .25 *( R->get_singlePos(particleNumber) - R->get_singlePos(i) )
                                    /((1+b*rr)*(1+b*rr)*rr);
        }


        for (int i=nParticleshalf;i<nParticles;i++)
            {
                    rr=R->get_rr(i-1,particleNumber);
                    result += .5 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                            /((1+b*rr)*(1+b*rr)*rr);
            }
    }

    else
    {
        for (int i=0;i<nParticleshalf;i++)
            {
                    rr=R->get_rr(particleNumber-1,i);
                    result += .5 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                            /((1+b*rr)*(1+b*rr)*rr);
            }
        for (int i=nParticleshalf;i<particleNumber;i++)
            {
            rr=R->get_rr(particleNumber-1,i);
            result += 0.25 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                    /((1+b*rr)*(1+b*rr)*rr);
            }

        for (int i=particleNumber+1;i<nParticles;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += .25 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                    /((1+b*rr)*(1+b*rr)*rr);
        }
    }

    return result;
}






vec TrialFct_molecule::GradSlater(int particleNumber, positions * R)
{
    vec result = zeros(ndim);
    result += gradhydrogen(particleNumber,0,R)/hydrogen(particleNumber,0,R);


    return result;
}




double TrialFct_molecule::ParamDerivativeOverFct(positions * R,int parameterNumber)
{
    double result =0.0;

    // if parameter is alpha, calculate only Slaterdet part
    if (parameterNumber==0)
    {
        for(int i=0;i<nParticles;i++)
        {
            result += dhydrogenda(i,0,R)/hydrogen(i,0,R);
        }

    }

    // if parameter is beta, calculate only Jastrow factor part
    else if(parameterNumber ==1)
    {
        double rr;
        for (int i=0;i<nParticleshalf;i++)
        {
            for(int j=0;j<i;j++)
            {
                rr = R->get_rr(i-1,j);
                result -=   0.25*rr/(  (1+funcParameters[parameterNumber]*rr)
                                 *(1+funcParameters[parameterNumber]*rr) ) ;

                rr = R->get_rr(i+nParticleshalf-1,j+nParticleshalf);
                result -=   0.25*rr/(  (1+funcParameters[parameterNumber]*rr)
                                 *(1+funcParameters[parameterNumber]*rr) ) ;
            }

            for (int j=nParticleshalf;j<nParticles;j++)
            {
                rr= R-get_rr(j-1,i);
                result -= 0.5*rr/(  (1+funcParameters[parameterNumber]*rr)
                                    *(1+funcParameters[parameterNumber]*rr) ) ;
            }
        }




    }



    // if parameter is neither alpha, beta or R_0 there is something wrong.
    // might be the distance between the nuclei in molecule calculations though.
    else
    {
        cout <<"ParamDerivativeOvecFct::Unknown parameter number"<<endl;
    }

    return result;
}
