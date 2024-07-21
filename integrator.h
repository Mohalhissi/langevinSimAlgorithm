void integrate(double &x, double &fx, double &en, \
               double &vx, double &randz, double &randth, \
               double &gamma, double &sigma, double &etot, \
               double &dt, double &mass, double &Eb){
    double A,fxc,xx,vvx;

    fxc=0;

    A=0.5*pow(dt,2)*(fx- gamma * vx)+ \
            sigma * pow(dt,(3.0/2))*(0.5*randz + \
                                     (0.5/sqrt(3.0))*randth);


    xx=x+(dt * vx)+A ;

    fxc= ((-1) * Eb/4.0)*pow(xx,3) + Eb * xx ;

    fxc=fxc/(mass) ;

    vvx=vx+(0.5*dt * (fxc+ fx))- dt * gamma * vx \
            + sigma*sqrt(dt) * randz - gamma * A ;

    etot=en+(0.5*pow(vvx,2)) ;

    x=xx ;
    vx=vvx ;
}

