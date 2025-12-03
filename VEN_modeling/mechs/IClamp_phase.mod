: phase current clamp

NEURON {
   POINT_PROCESS IClamp_phase
   RANGE del, dur, amp, A, delPhase, i
   ELECTRODE_CURRENT i
}

UNITS { (nA) = (nanoamp) }

PARAMETER {
   del	(ms)
   dur	(ms)  < 0, 1e9 >
   amp	(nA)
   A    (nA)  < 0, 10  >
   delPhase 
}

ASSIGNED { i (nA) }

INITIAL { i = 0 }

BREAKPOINT {
   at_time(del)
   at_time(del+dur)
   if ( t < del + dur && t > del ) {       
       i = amp+A*sin(6*3.1415926*t+delPhase)
   } else {
       i = 0
   }
}




