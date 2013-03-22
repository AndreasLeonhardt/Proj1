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


        // jastrow factor

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

        //result *= exp(jastrow); disable Jastrow factor for a moment

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
                    inverseSlaterUp(k,j)/=ratio;
                }
                else
                {   sum=0.0;
                    for (int l=0;l<nParticleshalf;l++)
                    {
                        sum+= hydrogen(particleNumber,l,Rnew)*inverseSlaterUp(l,j);
                    }
                    inverseSlaterUp(k,j) += -inverseSlaterUp(k,particleNumber)*sum/ratio;
                }
            }
        }
    }
    // if it has spind down, just update spin down
    else
    {
        for (int k=0;k<nParticleshalf;k++)
        {
            for (int j=0;j<nParticleshalf;j++)
            {
                if(particleNumber==j+nParticleshalf)
                {
                    inverseSlaterDown(k,j)/=ratio;
                }
                else
                {   sum=0.0;
                    for (int l=0;l<nParticleshalf;l++)
                    {
                        sum+= hydrogen(particleNumber,l,Rnew)*inverseSlaterDown(l,j);
                    }
                    inverseSlaterDown(k,j)+=-inverseSlaterDown(k,particleNumber-nParticleshalf)/ratio*sum;
                }
            }
        }
    }

}


//================================================================================================
// (laplace |D|)/|D|
double TrialFct_analytical::getDivGradOverFct(int particleNumber, positions *R)
{
    double result =0.0;


    // gradient
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

// with Jastrow factor: gradient of Jastrow factor and product of the gradients.



    return result;
}








//===========================================================================================
// returns the quantum force. using closed form expressions
// and the inverse slater determinant matrix. Be carefull, that this matrix is
// set to the right value.
vec TrialFct_analytical::quantumForce(int particleNumber, positions *R)
{
    vec result = zeros(ndim);
    //double beta = funcParameters[1];


    // particle has spin up
    if (particleNumber<nParticleshalf)
    {   // loop over spin up
        for (int j=0;j<nParticleshalf;j++)
        {
            // add slater determinant up part
            result += gradhydrogen(particleNumber,j,R)*inverseSlaterUp(j,particleNumber);
            // add jastrow factor part
            if(j!=particleNumber)
            {
               // results += rhat/(1+beta*r)/(1+beta*r)*0.5; skip Jastrow factor
            }
        }
        // jastrow factor for spin down particles
        for (int j=nParticleshalf;j<nParticles;j++)
        {
            //results += rhat/(1+beta*r)/(1+beta*r)*.25; skip jastrow factor
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
            // add Jastrow factor part with oppositen spin
          //  results += rhat/(1+beta*r)/(1+beta*r)*.25;
        }
        // Jastrow factor for particles with the same spin
        for(int j=nParticleshalf;j<nParticles;j++)
        {
            if(j!=particleNumber)
            {
            //    results += rhat/(1+beta*r)/(1+beta*r)*0.5;
            }
        }
    }

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


//===========================================================================================================
// hydrogen like wave functions, with parameter alpha
// hard coded for the first levels up to n=2, l=1, m_l = -1,0,1
// not including spin degeneracy.
double TrialFct_analytical::hydrogen(int particleNumber, int orbital, positions * R)
{
    double r=R->get_r(particleNumber);

    double result;

    if (orbital==0)
    {
        result = exp(-funcParameters[0]*r);
    }

    else if(orbital==1)
    {
        result = ( 1.0 - 0.5* funcParameters[0] * r )
                 * exp( -0.5* funcParameters[0] * r );
    }
    else if(orbital==2 || orbital==3 || orbital==4)
    {
        result = funcParameters[0]*R->get_singlePos(particleNumber)(orbital-2)
                *exp(-0.5*funcParameters[0]*r);
    }
    else
    {
        cout << "single particle wave function undefined"<<endl;
    }

    return result;
}

// gradient of hydrogen like wave functions
// hard coded as above
vec TrialFct_analytical::gradhydrogen(int particleNumber, int orbital, positions *R)
{
    vec result =zeros(ndim);
    double a = funcParameters[0];
    double r = R->get_r(particleNumber);

    if (orbital==0)
    {
        result = R->get_singlePos(particleNumber);
        result *=  -a/r *exp(-a*r);
    }

    else if(orbital==1)
    {
        result = R->get_singlePos(particleNumber);
        result *= a/(4*r)*(a*r-4)*exp(-a*r/2);
    }

    else if(orbital==2 || orbital==3 || orbital==4)
    {
        result = R->get_singlePos(particleNumber)*a*R->get_singlePos(particleNumber)(orbital);
        result(orbital)+= -2*r;
        result *= -a/(2*r)*exp(-a*r/2);
    }


    return result;
}


//laplace on hydrogen like wave functions
// closed form expressions, hard coded
double TrialFct_analytical::divgradhydrogen(int particleNumber, int orbital, positions* R)
{

    double result;
    double a = funcParameters[0];
    double r = R->get_r(particleNumber);

    if (orbital==0)
    {

        result = a/r*(a*r-2) * exp(-a*r);
    }

    else if(orbital==1)
    {
        result = -a/(8*r)*(a*a*r*r-10*a*r+18)*exp(-a*r/2);
    }

    else if(orbital==2 || orbital==3 || orbital==4)
    {
        vec RR = R->get_singlePos(particleNumber);
        RR *= a*a/(4*r)*(a*r-8)*exp(-a*r/2);

        result = RR(orbital);
    }

    return result;

}
