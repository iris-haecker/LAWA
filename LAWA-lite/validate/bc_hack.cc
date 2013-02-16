int
main()
{


   // trial space
   IntervalBasis   uTime(timeD, timeD_, timeJ0);
   uTime.enforceBoundaryCondition<NoBC>();
   
   // test space(s)
   IntervalBasis   vTime0(timeD, timeD_, timeJ0);
   IntervalBasis   vTime1(timeD, timeD_, timeJ0);
   
   vTime0.enforceBoundaryCondition<DirichletBC>();
   vTime1.enforceBoundaryCondition<NoBC>();
   
   
   
   
   Integral<Gauss,IntervalBasis,IntervalBasis>  integral0(uTime, vTime0);
   Integral<Gauss,IntervalBasis,IntervalBasis>  integral1(uTime, vTime1);
   
   
   if (kv<kvm(_jv)) {
   //  left of middle
       value = -integral1(ju, ku, uType, 0,
                          jv, kv, vType, 1)
               +integral1(ju, ku, uType, 0,
                          jv, kv, vType, 0);
   } else {
   //  right of middle
       value = -integral0(ju, ku, uType, 0,
                          jv, kv, vType, 1)
               +integral0(ju, ku, uType, 0,
                          jv, kv, vType, 0);
   }




}