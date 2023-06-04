README
In the “main” function there are three functions:
1.	The function “keySboxRecoveryAttackZD”, which performs the attack described in Section 4, as follows:
	a.	First, it generates all the needed variables and all the good pairs for the first stage of the attack.
	b.	The function “SboxRecoveryAttack0” recovers two S-boxes, as described in Section 4.2.
	c.	The function “fndAddS” recovers two additional S-boxes (the second stage), as described in Section 4.3.
	d.	The function “createThirdPCSet” creates good pairs for the third stage, and then the third and fourth stages (Sections 4.4 and 4.5) are performed in the function “keySboxRecoveryAttackZD” itself, for each remaining subkey.
	e.	Some notes:
		i.	The number of good pairs is defined at the beginning of the file under the name NGOODPAIRS (256 by default).
		ii.	The function “buildSbox” performs the strategy described in Section 4.1.
2.	The function “CollisionInMDGost”, which creates collisions in the GOST-based Merkle-Damgard hash function, for given IV and S-boxes.
3.	Two examples of finding collision, using specific Sboxes and IV=0.
