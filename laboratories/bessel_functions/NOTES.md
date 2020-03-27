
I tried to have a modern fancy CF implementation and Temme sum thing for sf_bessel.tcc and sf_mod_bessel.tcc.

I failed.

Actually, the modified Bessels worked just swimmingly. The usual J, N failed.

The Temme nk sum might still work for N.

In any case I need to roll back to bessels that work in the main library.
I will carry on, in a new branch, this bessel work.

The deal with the NaNs at z=2 is the cutover at that point to another algo - the pq CF.
The deal was that q was negative and blowing a sqrt.
Look for FIXME in sf_bessel ~711 after __cyl_hankel_1_ratio_j_frac.

Also, there is some deal with a sign for a CF.
