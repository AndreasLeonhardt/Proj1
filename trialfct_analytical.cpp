#include "trialfct_analytical.h"

TrialFct_analytical::TrialFct_analytical()
{
}
TrialFct_analytical::TrialFct_analytical(Config * parameters) : function(parameters)
{
}






// calculate the value of the trial function the position given in the variable position.
// Slater determinant of hydrogen wave functions and Jastrow factor.
double TrialFct_analytical::getValue(positions * R)
{

        double result = 0.0;

        // create Slater matrix
        mat Slaterup  = zeros(nParticleshalf,nParticleshalf);
        mat Slaterdown  = zeros(nParticleshalf,nParticleshalf);



        for (int i=0;i<nParticleshalf;i++)
        {
            for (int j=0;j<nParticleshalf;j++)
            {
                Slaterup(i,j)   = hydrogen(i,j,R);
                Slaterdown(i,j) = hydrogen(i+nParticleshalf,j,R);
            }
        }
        result = det(Slaterup)*det(Slaterdown);


       //  jastrow factor

        double jastrow = 0.0;
        double dist;

        for (int i=0;i<nParticleshalf;i++)
        {   // equal spin, factor 1/2
            for (int j=0;j<i;j++)
            {
                dist = R->get_rr(i-1,j);
                jastrow += dist/(2+2*funcParameters[1]*dist);

                dist = R->get_rr(nParticleshalf+i-1,nParticleshalf+j);
                jastrow += dist/(2+2*funcParameters[1]*dist);
            }

            // opposite spin, factor 1/4
            for (int j=nParticleshalf;j<nParticles;j++)
            {
                dist = R->get_rr(j-1,i);
                jastrow += dist/(4+4*funcParameters[1]*dist);
            }

        }

        result *= exp(jastrow);

        return result;
}



void TrialFct_analytical::setSlaterinv(positions * R)
{

    // create Slater matrix
    mat Slaterup  = zeros(nParticleshalf,nParticleshalf);
    mat Slaterdown  = zeros(nParticleshalf,nParticleshalf);



    for (int i=0;i<nParticleshalf;i++)
    {
        for (int j=0;j<nParticleshalf;j++)
        {
            Slaterup(i,j)   = hydrogen(i,               j,R);
            Slaterdown(i,j) = hydrogen(i+nParticleshalf,j,R);
        }
    }

    // calculate the inverse of the Slatermatrix for spin up and down
    // respectively. This is used for later updates of the position.
    inverseSlaterDown = inv(Slaterdown);
    inverseSlaterUp   = inv(Slaterup);

}


// after moving one particle, this updates the inverse of the Slater matrix,
// however, there are some issues, since the updatet matrix differs from the one
// set up from the new position with teh above function setSlaterinv.
void TrialFct_analytical::updateSlaterinv(int particleNumber, positions* Rnew, double ratio)
{
    mat newslatermat = zeros(nParticleshalf,nParticleshalf);
    double sum;
    // if particle has spin up, just update SlaterUp
    if(particleNumber<nParticleshalf)
    {
        for (int k=0;k<nParticleshalf;k++)
        {
            for (int j=0;j<nParticleshalf;j++)
            {
                if(particleNumber==j)
                {
                    newslatermat(k,j)=inverseSlaterUp(k,j)/ratio;
                }
                else
                {   sum=0.0;
                    for (int l=0;l<nParticleshalf;l++)
                    {
                        sum+= hydrogen(particleNumber,l,Rnew)*inverseSlaterUp(l,j);
                    }
                    newslatermat(k,j) =inverseSlaterUp(k,j) -inverseSlaterUp(k,particleNumber)*sum/ratio;
                }
            }
        }
        inverseSlaterUp = newslatermat;
    }
    // if it has spin down, just update spin down
    else
    {
        for (int k=0;k<nParticleshalf;k++)
        {
            for (int j=0;j<nParticleshalf;j++)
            {
                if(particleNumber==j+nParticleshalf)
                {
                    newslatermat(k,j)=inverseSlaterDown(k,j)/ratio;
                }
                else
                {   sum=0.0;
                    for (int l=0;l<nParticleshalf;l++)
                    {
                        sum+= hydrogen(particleNumber,l,Rnew)*inverseSlaterDown(l,j);
                    }
                    newslatermat(k,j)=inverseSlaterDown(k,j)-inverseSlaterDown(k,particleNumber-nParticleshalf)/ratio*sum;
                }
            }
        }
        inverseSlaterDown = newslatermat;
    }
}


//================================================================================================
// (laplace \Psi)/\Psi
// This includes several parts:
//   (laplace |D|) / |D|
// + div((grad Jastrow)/Jastrow)
// + ((grad Jastrow)/Jastrow)^2
// + 2* (grad |D|)/|D| * (grad Jastrow)/Jastrow

double TrialFct_analytical::getDivGradOverFct(int particleNumber, positions *R)
{
    double result =0.0;


    // div grad Slaterdet -------------------------------------------------------
    if(particleNumber<nParticleshalf)
    {
        for(int j=0;j<nParticleshalf;j++)
        {
            result += divgradhydrogen(particleNumber,j,R)*inverseSlaterUp(j,particleNumber);
            // calculate gradient here, if Jastrow
        }
    }
    else
    {
        for(int j=0;j<nParticleshalf;j++)
        {
            result += divgradhydrogen(particleNumber,j,R)*inverseSlaterDown(j,particleNumber-nParticleshalf);
            // calculate gradient here, if Jastrow included.
        }
    }

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
            result += 1 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }

        for (i=particleNumber+1;i<nParticleshalf;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += 1 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }

        for (i=nParticleshalf;i<nParticles;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += .5 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }
    }

    else
    {
        int i;
        for (int i=0;i<nParticleshalf;i++)
        {
            rr=R->get_rr(particleNumber-1,i);
            result += .5 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }
        for (i=nParticleshalf;i<particleNumber;i++)
        {
            rr=R->get_rr(particleNumber-1,i);
            result += 1 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
        }

        for (i=particleNumber+1;i<nParticles;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += 1 / ((1+b*rr)*(1+b*rr)) * (  1/rr - b/(1+b*rr)  );
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
vec TrialFct_analytical::quantumForce(int particleNumber, positions *R)
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
double TrialFct_analytical::SlaterRatio(int particleNumber ,positions * Rold,positions * Rnew)
{    double result = 0;

    // calculate spin up part, if only spin up particle moved
    if(particleNumber<nParticleshalf)
    {

        for (int j=0;j<nParticleshalf;j++)
        {
            result += hydrogen(particleNumber, j, Rnew)*inverseSlaterUp(j,particleNumber);
        }
    }
    // calculate spin down ration, if only spin down particle moved
    else
    {
        for (int j=0;j<nParticleshalf;j++)
        {
            result += hydrogen(particleNumber, j, Rnew)
                    *inverseSlaterDown(j,particleNumber-nParticleshalf);
        }
    }

    return result;
}


double TrialFct_analytical::JastrowRatio(int particleNumber, positions * Rold, positions * Rnew)
{
    // jastrow factor new
    double jastrow = 0.0;
    double dist;



    if (particleNumber<nParticleshalf)
    {
        for (int i=0;i<particleNumber;i++)
        {
            dist = Rnew->get_rr(particleNumber-1,i);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = Rold->get_rr(particleNumber-1,i);
            jastrow -= dist/(2+2*funcParameters[1]*dist);
        }
        for (int i=particleNumber+1;i<nParticleshalf;i++)
        {
            dist = Rnew->get_rr(i-1,particleNumber);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = Rold->get_rr(i-1,particleNumber);
            jastrow -= dist/(2+2*funcParameters[1]*dist);
        }
        for (int i=nParticleshalf;i<nParticles;i++)
        {
            dist = Rnew->get_rr(i-1,particleNumber);
            jastrow += dist/(4+4*funcParameters[1]*dist);

            dist = Rold->get_rr(i-1,particleNumber);
            jastrow -= dist/(4+4*funcParameters[1]*dist);
        }
    }
    else
    {
        for (int i=0;i<nParticleshalf;i++)
        {
            dist = Rnew->get_rr(particleNumber-1,i);
            jastrow += dist/(4+4*funcParameters[1]*dist);

            dist = Rold->get_rr(particleNumber-1,i);
            jastrow -= dist/(4+4*funcParameters[1]*dist);
        }
        for (int i=nParticleshalf;i<particleNumber;i++)
        {
            dist = Rnew->get_rr(particleNumber-1,i);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = Rold->get_rr(particleNumber-1,i);
            jastrow -= dist/(2+2*funcParameters[1]*dist);
        }
        for (int i=particleNumber+1;i<nParticles;i++)
        {
            dist = Rnew->get_rr(i-1,particleNumber);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = Rold->get_rr(i-1,particleNumber);
            jastrow -= dist/(2+2*funcParameters[1]*dist);
        }
    }



    return exp(jastrow);
}

// calculates the divergence of the Jastrow factor divided by the jastrow factor
// That is sum_{i!=k} \frac{a}{ (1+\beta*r_{ki})^2 }  \frac{ \vec{r}_{ki} }{ r_{ki} }
// where a is 0.5 for spin_i = spin_k and
//            0.25 for spin_i!=spin_k
vec TrialFct_analytical::GradJastrow(int particleNumber, positions * R)
{
    vec result = zeros(ndim);
    double b = funcParameters[1];
    double rr;
    if(particleNumber<nParticleshalf)
    {
        for (int i=0;i<particleNumber;i++)
            {
            rr=R->get_rr(particleNumber-1,i);
            result += 0.5 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                    /((1+b*rr)*(1+b*rr)*rr);
            }

        for (int i=particleNumber+1;i<nParticleshalf;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += .5 *( R->get_singlePos(particleNumber) - R->get_singlePos(i) )
                                    /((1+b*rr)*(1+b*rr)*rr);
        }


        for (int i=nParticleshalf;i<nParticles;i++)
            {
                    rr=R->get_rr(i-1,particleNumber);
                    result += .25 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                            /((1+b*rr)*(1+b*rr)*rr);
            }
    }

    else
    {
        for (int i=0;i<nParticleshalf;i++)
            {
                    rr=R->get_rr(particleNumber-1,i);
                    result += .25 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                            /((1+b*rr)*(1+b*rr)*rr);
            }
        for (int i=nParticleshalf;i<particleNumber;i++)
            {
            rr=R->get_rr(particleNumber-1,i);
            result += 0.5 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                    /((1+b*rr)*(1+b*rr)*rr);
            }

        for (int i=particleNumber+1;i<nParticles;i++)
        {
            rr=R->get_rr(i-1,particleNumber);
            result += .5 *(R->get_singlePos(particleNumber)-R->get_singlePos(i))
                                    /((1+b*rr)*(1+b*rr)*rr);
        }
    }

    return result;
}






vec TrialFct_analytical::GradSlater(int particleNumber, positions * R)
{
    vec result = zeros(ndim);

    // particle has spin up
    if (particleNumber<nParticleshalf)
    {   // loop over spin up
        for (int j=0;j<nParticleshalf;j++)
        {
            // add slater determinant up part
            result += gradhydrogen(particleNumber,j,R)*inverseSlaterUp(j,particleNumber);
        }

    }

    // particle has spin down
    else
    {
        for(int j=0;j<nParticleshalf;j++)
        {
            // add slater determinant down part
            result += gradhydrogen(particleNumber,j,R)
                       *inverseSlaterDown(j,particleNumber-nParticleshalf);
        }
    }

    return result;
}




double TrialFct_analytical::ParamDerivativeOverFct(positions * R,int parameterNumber)
{
    double result =0.0;

    // if parameter is alpha, calculate only Slaterdet part
    if (parameterNumber==0)
    {
        for (int i=0;i<nParticleshalf;i++)
        {
            for(int j=0;j<nParticleshalf;j++)
            {
                result+=inverseSlaterUp(i,j)*dhydrogenda(j,i,R);
                result+=inverseSlaterDown(i,j)*dhydrogenda(j+nParticleshalf,i,R);
            }
        }

    }

    // if parameter is beta, calculate only Jastrow factor part
    else if(parameterNumber ==1)
    {
        double rr;
        for (int i=0;i<nParticles;i++)
        {
            for(int j=0;j<i;j++)
            {
                rr = R->get_rr(i-1,j);
                result -=   rr/(  (1+funcParameters[parameterNumber]*rr)
                                 *(1+funcParameters[parameterNumber]*rr) ) ;
            }
        }
    }

    // if parameter is neither alpha or beta, there is something wrong.
    // might be the distance between the nuclei in molecule calculations though.
    else
    {
        cout <<"ParamDerivativeOvecFct::Unknown parameter number"<<endl;
    }

    return result;
}
