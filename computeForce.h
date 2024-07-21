void force(double &fx, double &en, double &x, \
           double &mass, double &Eb){

    en=0 ; fx=0;

    fx= ((-1)*Eb/4.0)*pow(x,3) + Eb*x;

    fx=fx/(mass) ;

    en= (Eb/16)*pow(x,4)- (Eb/2)*pow(x,2) + Eb ;


}
