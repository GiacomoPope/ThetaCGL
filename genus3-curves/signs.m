/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */


function booltoGF2(bool)
  if bool then
  return GF(2)!1;
  else return GF(2)!0;
end if;
end function;

function chartoInt(cha)
   seq:=Eltseq(cha);
   retpre:=&+[2^(6-i) * Integers()!(GF(2)!seq[i]): i in [1..6]];
   return ((retpre -1) mod 64)+1;
end function;

function preloop()
zer:=ZeroMatrix(GF(2), 3,3);
id:=IdentityMatrix(GF(2),3);
J:=BlockMatrix(2,2, [zer,id, id,zer]);
J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
chars:=EvenThetaCharacteristics(3);

V:=VectorSpace(GF(2), 21);
cs:=[[GF(2)!1: i in [1..36]]];
vec:=[Matrix([[GF(2)!c: c in cha[1] cat cha[2]]]): cha in chars]; 
for v in Basis(V) do
	T:=UpperTriangularMatrix(Eltseq(v));
new:= [(ve*T*Transpose(ve))[1,1]: ve in vec];
        Append(~cs, new);
end for;
A:=Matrix(cs);
Ech:=EchelonForm(A);
pivots:=[Min([i: i in [1..36]| Ech[j,i] eq GF(2)!1]): j in [1..21]];
pivotscompl:=[i: i in [1..36]| not (i in pivots)];
W:=VectorSpace(GF(2), 6);
T:=[ W!(chars[i][1] cat chars[i][2]) : i in pivots];
Tcompl:=[ W!(chars[i][1] cat chars[i][2]) : i in pivotscompl];
bs:=&cat[[[b1,b2]: b1 in W| (Matrix(b1)*J*Transpose(Matrix(b2)))[1,1] eq 0 and Position(Eltseq(b1), GF(2)!1) lt Position(Eltseq(b2), GF(2)!1) and Eltseq(b1)[Position(Eltseq(b2), GF(2)!1)] eq GF(2)!0 and b1 ne W!0 and b2 ne W!0  and  (Matrix(b1)*J1*Transpose(Matrix(b1)))[1,1] eq (Matrix(b2)*J1*Transpose(Matrix(b2)))[1,1]   ]: b2 in W];
N:=#bs;
cs := [[chars[i] : i in [1..36]|  (Matrix(b[1])*J*Transpose(Matrix(vec[i])))[1,1] eq (Matrix(b[1])*J1*Transpose(Matrix(b[1])))[1,1] and  (Matrix(b[2])*J*Transpose(Matrix(vec[i])))[1,1] eq (Matrix(b[2])*J1*Transpose(Matrix(b[2])))[1,1] ]: b in bs  ];
rep := [[c  : c in cs[i] | (c[1] cat c[2])[Position(Eltseq(bs[i][1]), GF(2)!1)] eq 0 and (c[1]cat c[2])[Position(Eltseq(bs[i][2]), GF(2)!1)] eq 0] : i in [1..#bs]];
cosets := [[ [W!(rep[i][j][1] cat rep[i][j][2]), (W!(rep[i][j][1] cat rep[i][j][2])+bs[i][1]),  (W!(rep[i][j][1] cat rep[i][j][2])+bs[i][2]) , (W!(rep[i][j][1] cat rep[i][j][2])+bs[i][1]+bs[i][2]) ]: j in [1..3]   ]  : i in [1..N]];
S1:= [i  : i in [1..N] | &and[cosets[i][1][j] in T: j in [1..4] ]];
S2:= [i  : i in [1..N] | &and[cosets[i][2][j] in T: j in [1..4] ]];
S3:= [i  : i in [1..N] | &and[cosets[i][3][j] in T: j in [1..4] ]];
Si:=S1 cat S2 cat S3;

A1:=[Matrix([ [booltoGF2( &or[cosets[i][2][j] eq c: j in [1..4]])  : c in Tcompl], [booltoGF2(&or[cosets[i][3][j] eq c: j in [1..4]])       : c in Tcompl]]):  i in S1];
A2:=[Matrix([ [booltoGF2( &or[cosets[i][1][j] eq c: j in [1..4]])  : c in Tcompl], [booltoGF2(&or[cosets[i][3][j] eq c: j in [1..4]])       : c in Tcompl]]):  i in S2];
A3:=[Matrix([ [booltoGF2( &or[cosets[i][1][j] eq c: j in [1..4]])  : c in Tcompl], [booltoGF2(&or[cosets[i][2][j] eq c: j in [1..4]])       : c in Tcompl]]):  i in S3];
Ai:=A1 cat A2 cat A3;
As:=A1[1];
i:=2;
list:=[Si[1]];
listcoind:=[1];
while Rank(As) lt 15 do
       VJ:=VerticalJoin(As, Ai[i]);
       if Rank(VJ) gt Rank(As) then
 	 As:=VJ;
	 Append(~list, Si[i]);
	 if i lt #S1 then
	     Append(~listcoind, 1);
	 elif i lt #S1+#S2 then
	     Append(~listcoind, 2);
	 else
	     Append(~listcoind, 3);
         end if;
       end if;
       i:=i+1;
end while;
usedbs:=[bs[i]: i in list];
usedrep:=[rep[i]: i in list];
return As, usedbs, usedrep, listcoind, Tcompl;
end function;



intrinsic correct_signs(thetas::SeqEnum)-> SeqEnum
{corrects the signs of thetas}
X, bs, rep, coind, Tcompl:=preloop();
ZZ:=Integers();
W:=Parent(bs[1][1]);
W1:=RSpace(ZZ,6);
zer:=ZeroMatrix(ZZ, 3,3);
id:=IdentityMatrix(ZZ,3);
J:=BlockMatrix(2,2, [zer,id, id,zer]);
J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
N:=#bs;
cosets := [[ [W1!(rep[i][j][1] cat rep[i][j][2]), (W1!(rep[i][j][1] cat rep[i][j][2])+W1!bs[i][1]),  (W1!(rep[i][j][1] cat rep[i][j][2])+W1!bs[i][2]) , (W1!(rep[i][j][1] cat rep[i][j][2])-W1!bs[i][1]-W1!bs[i][2]) ]: j in [1..3]   ]  : i in [1..N]];
cosets4:=[[[ [cosets[i][j][nu], cosets[i][j][nu]- W1!bs[i][1], cosets[i][j][nu]- W1!bs[i][2], cosets[i][j][nu]- W1!bs[i][1]- W1!bs[i][2] ]: nu in [1..4]  ]: j in [1..3]]: i in [1..N]];
carry:=[[[[W1![Eltseq(cosets4[i][j][xi][mu])[nu] div 2: nu in [1..6]]:mu in [1..4]]: xi in [1..4]]: j in [1..3]]     : i in [1..N]];
vec:=[];
for i in [1..N] do
        // Term of the theta relation where the sign is fixed.
        coindi:=coind[i];
        // Other two terms
        compl:=[j: j in [1..3]| j ne coindi];
        coeff:=[&+[  (-1)^ &+[(Matrix(cosets4[i][j][xi][mu])*J1*Transpose(Matrix(carry[i][j][xi][mu]) )) [1,1]: xi in [1..4]]: mu in [1..4]  ]: j in [1..3]];
        if rep[i][1]  eq [[0,0,0],[0,0,0]] then
		coeff[1]-:=8;
	end if;
	sigsposs:=[[coeff[j]: nu in [1..4]]: j in [1..3]  ];
	/* nu=1 means no sign switch
	   nu=2 means 1st sign switches
	   nu=3 means 2nd sign switches
	   nu=4 means both signs switch
	*/
	sigsposs[compl[1]][2] *:= (-1);
	sigsposs[compl[2]][3] *:= (-1);
	sigsposs[compl[1]][4] *:= (-1);
	sigsposs[compl[2]][4] *:= (-1);
	charsInt:=[[chartoInt(cosets[i][j][xi]   )  : xi in [1..4]]: j in [1..3]];
	thetarel:=[Abs(&+[sigsposs[j][nu]*  &*[thetas[charsInt[j][xi] ]: xi in [1..4   ]]: j in [1..3]])   : nu in [1..4]];
	min, ind:=Min(thetarel);
        vec cat:= Intseq(ind-1, 2,2);	
end for;
v:=Vector(GF(2), vec);
sol:=Eltseq(Solution(Transpose(X), v));
for j in [1..#Tcompl ] do
      indi:=chartoInt(Tcompl[j]);
      thetas[indi] *:= (-1)^(ZZ!sol[j]);
end for;
return thetas;
end intrinsic;
