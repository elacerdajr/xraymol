#!/usr/bin/perl

use Math::Trig;

#Script que le as coordenadas cartesianas de uma molecula e fornece sua topologia.
#autor: Evanildo Laceda Jr.
#Intituto de Fisica, USP, Brasil. 

###################################
#Ultima modificacao:  27/10/2010 #
##################################
@comand = @ARGV;

$key =grep(/-key/,@comand);
$itp =grep(/-itp/,@comand);
$help =grep(/-h/,@comand);


if ($help==1){

&intro();
exit;

}

if (($key == 0) && ($itp == 0)){chomp($arquivo1 = "$ARGV[0]");}
else {chomp($arquivo1 = "$ARGV[1]");}



if ($ARGV[2] != Null){$maxbondlen = $ARGV[2];}
else {$maxbondlen = 1.56;}


@arq1 = split(/\./,$arquivo1);
$ext=$arq1[1];

open(in,"<$arquivo1") || die "$!: $arquivo1";

@linhas = <in>;

if ($ext eq "xyz") {

chomp($N=$linhas[0]);
$N =~ s/ +//g;

@caixa =split(/=/,$linhas[1]);
chomp($caixa[1]);
#print "$caixa[1]\n";

if ($key == 1) { $arq_xray = $arq1[0] . ".key";}
if ($itp == 1) { $arq_xray = $arq1[0] . ".itp";}
else {$arq_xray = $arq1[0] . ".xray"}

open(out,">$arq_xray");



for ($i=1; $i<=$N; $i++) 
{

@atom = split(/ +/,$linhas[$i+1]);
chomp(@atom);


$sym[$i] = $atom[0];
$x[$i] = $atom[1]; 
$y[$i] = $atom[2]; 
$z[$i] = $atom[3];

#printf "%s     %6.3f %6.3f %6.3f\n",$sym[$i],$x[$i],$y[$i],$z[$i]; 

}


########################
#identificando ligacoes#
########################


$n_bon = 0;

 for ($i=1; $i<=$N; $i++) {
   
   for ($j=$N; $j>=$i+1; $j--) {

    $x_2 = ($x[$i]-$x[$j])*($x[$i]-$x[$j]);
    $y_2 = ($y[$i]-$y[$j])*($y[$i]-$y[$j]);
    $z_2 = ($z[$i]-$z[$j])*($z[$i]-$z[$j]);     

    $d = sqrt($x_2 + $y_2 + $z_2);

    if ($d < $maxbondlen) { 
      $bnd = sprintf("%4d  %4d      %5.3f\n",$i,$j,$d);
      $bnd_itp = sprintf("%4d  %4d      %7.4f\n",$i,$j,$d/10);
      push(@bond,$bnd);
      push(@bond_itp,$bnd_itp);
      $n_bon++;
      $k[$i]++;
      $connection[$i][$j] = "    $j";
      $connection[$j][$i] = "    $i";
      }
    else{
        $connection[$i][$j] = "";
        $connection[$j][$i] = "";
        }
    }
   
   

   }



#########################
# identificando angulos #
#########################

$fmt_ang= "%4d  %4d  %4d      %9.2f\n";

for ($i=0; $i<=$#bond; $i++) 
{ 
for ($j=$i+1; $j<=$#bond; $j++) 
{
@at_i= split(/ +/,$bond[$i]);
@at_j= split(/ +/,$bond[$j]);

if ($at_i[2] == $at_j[1])  {
$a = $at_i[1]; $b = $at_i[2]; $c=$at_j[2];
$theta = &ang($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c]); 
push(@angle, sprintf($fmt_ang,$a,$b,$c,$theta));
#push(@angle,sprintf("%4d  %4d  %4d\n", $at_i[0],$at_i[1],$at_j[1]));
} 

elsif ($at_i[1] == $at_j[1]) { 
$a = $at_j[2]; $b = $at_i[1]; $c=$at_i[2];
$theta = &ang($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c]); 
push(@angle, sprintf($fmt_ang,$a,$b,$c,$theta));
}
 
elsif ($at_i[1] == $at_j[2]) { 
$a = $at_i[2]; $b = $at_i[1]; $c=$at_j[1];
$theta = &ang($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c]); 
push(@angle, sprintf($fmt_ang,$a,$b,$c,$theta));
#push(@angle,sprintf("%4d  %4d  %4d\n",$at_i[1],$at_i[0],$at_j[0]));
} 


elsif  ($at_i[2] == $at_j[2]) {
$a = $at_i[1]; $b = $at_i[2]; $c=$at_j[1];
$theta = &ang($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c]); 
push(@angle, sprintf("%4d  %4d  %4d      %9.2f\n",$a,$b,$c,$theta));
#push(@angle,sprintf("%4d  %4d  %4d\n",$at_i[0],$at_i[1],$at_j[0]));
} 



}
}


@angle = sort(@angle);

#######################
#identificando torcoes#
#######################

$fmt_dhd= "%4d  %4d  %4d  %4d        %9.2f\n";

for ($i=0; $i<=$#angle; $i++) 
{ 
for ($j=$i+1; $j<=$#angle; $j++) 
{


@at_ang_i= split(/ +/,$angle[$i]);
@at_ang_j= split(/ +/,$angle[$j]);


 

    if (($at_ang_i[1] == $at_ang_j[2]) && ($at_ang_i[2] == $at_ang_j[3])) { 
      
      $a=$at_ang_j[1]; $b=$at_ang_j[2]; $c=$at_ang_j[3]; $d=$at_ang_i[3];
      $dhd = &torcao($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c],$x[$d],$y[$d],$z[$d]); 
      $tor=sprintf($fmt_dhd,$a,$b,$c,$d,$dhd);
      push(@torsion,$tor);
    } 
    elsif (($at_ang_i[2] == $at_ang_j[1]) && ($at_ang_i[3] == $at_ang_j[2])) {
      $a=$at_ang_i[1]; $b=$at_ang_i[2]; $c=$at_ang_i[3]; $d=$at_ang_j[3];
      $dhd = &torcao($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c],$x[$d],$y[$d],$z[$d]); 
      $tor=sprintf($fmt_dhd,$a,$b,$c,$d,$dhd);
      push(@torsion,$tor);
    } 
    elsif (($at_ang_i[1] == $at_ang_j[2]) && ($at_ang_i[2] == $at_ang_j[1])) {
      $a=$at_ang_i[3]; $b=$at_ang_i[2]; $c=$at_ang_i[1]; $d=$at_ang_j[3];
      $dhd = &torcao($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c],$x[$d],$y[$d],$z[$d]); 
      $tor=sprintf($fmt_dhd,$a,$b,$c,$d,$dhd);
      push(@torsion,$tor); 
#   $tor=sprintf("%4d  %4d  %4d  %4d\n",$at_ang_j[3],$at_ang_i[1],$at_ang_i[2],$at_ang_i[3]);
    } 
    elsif (($at_ang_i[2] == $at_ang_j[3]) && ($at_ang_i[3] == $at_ang_j[2])) {
      $a=$at_ang_i[1]; $b=$at_ang_i[2]; $c=$at_ang_i[3]; $d=$at_ang_j[1];
      $dhd = &torcao($x[$a],$y[$a],$z[$a],$x[$b],$y[$b],$z[$b],$x[$c],$y[$c],$z[$c],$x[$d],$y[$d],$z[$d]); 
      $tor=sprintf($fmt_dhd,$a,$b,$c,$d,$dhd);
      push(@torsion,$tor);
 #    $tor=sprintf("%4d  %4d  %4d  %4d\n",$at_ang_i[1],$at_ang_i[2],$at_ang_i[3],$at_ang_j[1]);
      } 

   }
}

@bond = sort(@bond); 
@bond_itp = sort(@bond_itp);
@angle = sort(@angle);
@torsion = sort(@torsion);

#######################
#Imprimindo resultados#
#######################

 $n_ang =$#angle +1;
 $n_tor = $#torsion +1;
 
&intro();

 print "\# atoms = $N\n";
 print "\# bonds = $n_bon\n";
 print "\# angles = $n_ang\n";
 print "\# torsions = $n_tor\n";
 print "\# bond criteria dij < $maxbondlen\n\n";
 



 
 print "input: " . $arquivo1 ."\n";
 print "output: $arq_xray\n\n";

# Arquivo .key do TINKER!

if ($key == 1){

print out "\n\n\# Output Control\nVERBOSE\n";
print out "\# atoms = $N\n";
print out "\# bonds = $n_bon\n";
print out "\# angles = $n_ang\n";
print out "\# torsions = $n_tor\n";
print out "\#bond criteria < $maxbondlen\n\n";

print out "\n\# Force Field Selection\n";
print out "PARAMETERS        none\n";

print out "CHARGETERM\n";
print out "VDWTERM\n";
print out "ANGLETERM\n";      
print out "BONDTERM\n";     
print out "TORSIONTERM\n";

print out "\n\# Van Der Waals Functional Form\n";
print out "EPSILONRULE        GEOMETRIC\n";
print out "RADIUSRULE         GEOMETRIC\n";
print out "RADIUSSIZE         DIAMETER\n";
print out "RADIUSTYPE         SIGMA\n";
print out "VDWTYPE            LENNARD-JONES\n";
print out "VDW-14-SCALE       2.0                \#OPLSAA\n";
print out "CHG-14-SCALE       2.0                \#OPLSAA\n";
print out "TORSIONUNIT        0.5                \#OPLSAA\n";

print out "\n\# Electrostatics Functional Form\n";
print out "DIELECTRIC        1.0\n";

print out "\n\# Dynamics\n";
print out "INTEGRATE          VERLET\n";
print out "TAU-TEMPERATURE   0.1\n";
print out "THERMOSTAT         BERENDSEN\n";

print out "\n\# Random Number\n";
print out "RANDOMSEED        12345789\n\n";


print out "\# Molecula 1 - natom = $N\n";
print out "\#ATOM    i    Sym    Descr        natm    mass   valence\n\n";
for ($i=1; $i<=$N; $i++) {
print out "atom     $i    $sym[$i]      \"$sym[$i]   mol1\"\n"; 
} 

print out "\n#####        Charges \n";
print out "\n#####      i    q\n\n";  
for ($i=1; $i<=$N; $i++) {
print out "charge    $i \n"; 
} 


print out "\n#####         Van der Waals Parameters \n";
print out "\n#####        i    sig     eps\n\n";  
for ($i=1; $i<=$N; $i++) {
print out "vdw    $i \n"; 
} 



print out "\n#####        Bonds \n";
print out "\n#####     i    j        r    \n\n";  

  for ($i=0; $i<=$#bond; $i++) { print out "bond   $bond[$i]";}

print out "\n#####        Angles \n";
print out "\n#####      i     j     k          th    \n\n";  

  for ($i=0; $i<=$#angle; $i++) { print out "angle   $angle[$i]";}
 
print out "\n#####        torsions \n";
print out "\n#####        i     j     k     l            dhd    \n\n";  

   for ($i=0; $i<=$#torsion; $i++) { print out "torsion   $torsion[$i]";}
} 

# .itp flie

if ($itp == 1){
print out "\; atoms = $N\n";
print out "\; bonds = $n_bon\n";
print out "\; angles = $n_ang\n";
print out "\; torsions = $n_tor\n";
print out "\; bond criteria < $maxbondlen\n\n";

print out "\n[ bonds ]\n\n"; 
  for ($i=0; $i<=$#bond_itp; $i++) { print out $bond_itp[$i];}

print out "\n[ angles ]\n\n"; 
  for ($i=0; $i<=$#angle; $i++) { print out $angle[$i];}
}

else { 
 print out "\# atoms = $N\n";
 print out "\# bonds = $n_bon\n";
 print out "\# angles = $n_ang\n";
 print out "\# torsions = $n_tor\n";
 print out "\#bond criteria < $maxbondlen\n\n";

print out "\n[ connections ]\n\n";
for ($i=1; $i<=$#connection; $i++) { 
    printf out "%3d",$i;
 for ($j=1; $j<=$#connection; $j++) { 
    print out $connection[$i][$j];
 }
 print out "\n";
}


print out "\n[ bonds ]\n\n"; 
  for ($i=0; $i<=$#bond; $i++) { print out $bond[$i];}

print out "\n[ angles ]\n\n"; 
  for ($i=0; $i<=$#angle; $i++) { print out $angle[$i];}
 
print out "\n[ torsions ]\n\n"; 
   for ($i=0; $i<=$#torsion; $i++) { print out $torsion[$i];}
}
 

 }

else {print "Erro: input invalido. Use a extensao *.xyz.\n"; exit;}


close(in);
close(out); 




# @nodes = map(m/(-key/g,$tor);
# print "\$tor = $tor";
# for ($i=0; $i<=$#nodes; $i++){
# print $nodes[$i] . "\n";
# }




#########################################################################################

sub intro {
 print "------------------------------------\n";
 print "X-raymol\n";
 print "created by Evanildo Lacerda Jr.\n";
 print "------------------------------------\n\n";
 
 print "command -> \$ xraymol [option] [file.xyz] \n\n";
 print "[options]  \n";
 print "none \t generate .xray file\n";
 print "-key \t generate .key file [tinker]\n";
 print "-itp \t generate .itp file [gromacs]\n";
 print "-h   \t help\n\n";
}

sub ang
{

$rad2deg=57.29577951;

$xi=$_[0];  $yi=$_[1]; $zi=$_[2];
$xj=$_[3];  $yj=$_[4]; $zj=$_[5];
$xk=$_[6];  $yk=$_[7]; $zk=$_[8];

# j -> i 

$xji=$xi - $xj;
$yji=$yi - $yj; 
$zji=$zi - $zj;

$dji = sqrt($xji*$xji + $yji*$yji + $zji*$zji);

#j -> k
$xjk=$xk - $xj;
$yjk=$yk - $yj; 
$zjk=$zk - $zj;

$djk = sqrt($xjk*$xjk + $yjk*$yjk + $zjk*$zjk);

$costheta = ($xji*$xjk + $yji*$yjk + $zji*$zjk)/($dji*$djk);

$theta = acos($costheta);

return $theta*$rad2deg;
}



sub torcao
{
$rad2deg=57.29577951;

$xi=$_[0];  $yi=$_[1]; $zi=$_[2];
$xj=$_[3];  $yj=$_[4]; $zj=$_[5];
$xk=$_[6];  $yk=$_[7]; $zk=$_[8];
$xl=$_[9];  $yl=$_[10]; $zl=$_[11];

# j -> i 

$xji=$xi - $xj;
$yji=$yi - $yj; 
$zji=$zi - $zj;

$dji = sqrt($xji*$xji + $yji*$yji + $zji*$zji);

#j -> k
$xjk=$xk - $xj;
$yjk=$yk - $yj; 
$zjk=$zk - $zj;

$djk = sqrt($xjk*$xjk + $yjk*$yjk + $zjk*$zjk);

$nxijk = $yjk*$zji - $zjk*$yji;
$nyijk = $zjk*$xji - $xjk*$zji;
$nzijk = $xjk*$yji - $yjk*$xji;

$nijk = sqrt($nxijk*$nxijk + $nyijk*$nyijk + $nzijk*$nzijk); 

#k -> l
$xkl=$xl - $xk;
$ykl=$yl - $yk; 
$zkl=$zl - $zk;

# k -> j

$xkj= -$xjk;
$ykj= -$yjk;
$zkj= -$zjk;


$nxjkl =  $ykl*$zkj - $zkl*$ykj;
$nyjkl =  $zkl*$xkj - $xkl*$zkj;
$nzjkl =  $xkl*$ykj - $ykl*$xkj;

$njkl = sqrt($nxjkl*$nxjkl + $nyjkl*$nyjkl + $nzjkl*$nzjkl); 

$escnn = ($nxijk*$nxjkl + $nyijk*$nyjkl + $nzijk*$nzjkl)/$nijk/$njkl;

$phi = acos($escnn)*$rad2deg;

$horario = &ehhorario($xi,$yi,$zi,$xj,$yj,$zj,$xk,$yk,$zk,$xl,$yl,$zl);

if ($horario) { $phi = 360-$phi;}
      
return $phi;      

}

sub ehhorario
{
$rad2deg=57.29577951;

my $m1,$m2;

$pa[1]=$_[0];  $pa[2]=$_[1]; $pa[3]=$_[2];
$pb[1]=$_[3];  $pb[2]=$_[4]; $pb[3]=$_[5];
$pc[1]=$_[6];  $pc[2]=$_[7]; $pc[3]=$_[8];
$pd[1]=$_[9];  $pd[2]=$_[10]; $pd[3]=$_[11];


# CEntrando em pb
for($m1=1; $m1<=3; $m1++){
$pa[$m1]= $pa[$m1] - $pb[$m1];
$pc[$m1]= $pc[$m1] - $pb[$m1];
$pd[$m1]= $pd[$m1] - $pb[$m1];
$pb[$m1]= $pb[$m1] - $pb[$m1];
}

#     Definindo o vetor do eixo de rotacao

$ux=$pc[1]; $uy=$pc[2]; $uz=$pc[3];

$dcb = sqrt($ux*$ux + $uy*$uy + $uz*$uz);


#     Angulo de rotacao em torno de u
      $gama = 0.01/$rad2deg;

#     Angulos de projecao            
      $hcb = sqrt($ux*$ux+$uy*$uy);
      $costheta = $ux/$hcb;
      $sintheta = $uy/$hcb;
      $cosphi   = $uz/$dcb;
      $sinphi   = $hcb/$dcb;

#     Matrizes de rotacao 
      $A[1][1] =  $costheta;
      $A[2][2] =  $costheta;
      $A[3][3] =  1;
      $A[1][2] =  $sintheta;
      $A[2][1] = -$sintheta;
      $A[1][3] =  0;
      $A[2][3] =  0;
      $A[3][1] =  0;
      $A[3][2] =  0;

      $B[1][1] =  $cosphi;
      $B[2][2] =  1;
      $B[3][3] =  $cosphi;
      $B[1][3] = -$sinphi;
      $B[3][1] =  $sinphi;
      $B[1][2] =  0;
      $B[2][1] =  0;
      $B[3][2] =  0;
      $B[2][3] =  0;

      $C[1][1] =  cos($gama);
      $C[2][2] =  cos($gama);
      $C[3][3] =  1;
      $C[1][2] =  sin($gama);
      $C[2][1] = - sin($gama);
      $C[1][3] =  0;
      $C[2][3] =  0;
      $C[3][1] =  0;
      $C[3][2] =  0;

#     Novas coordenadas

#      BA  = matmul(B,A)
    for($m1=1; $m1<=3; $m1++){
      for($m2=1; $m2<=3; $m2++){
        $BA[$m1][$m2]= $B[$m1][1]*$A[1][$m2]+$B[$m1][2]*$A[2][$m2]+$B[$m1][3]*$A[3][$m2];
       }
      }

#      CBA = matmul(C,BA)

    for($m1=1; $m1<=3; $m1++){
      for($m2=1; $m2<=3; $m2++){
        $CBA[$m1][$m2]= $C[$m1][1]*$BA[1][$m2]+$C[$m1][2]*$BA[2][$m2]+$C[$m1][3]*$BA[3][$m2];
       }
      }

#      panew = matmul(BA,pa)
#      pbnew = matmul(BA,pb)      
#      pcnew = matmul(BA,pc)
#      pdnew = matmul(BA,pd)
#      pdnewr = matmul(CBA,pd)
  
     for($m1=1; $m1<=3; $m1++){
        $panew[$m1]= $BA[$m1][1]*$pa[1]+$BA[$m1][2]*$pa[2]+$BA[$m1][3]*$pa[3];
        $pbnew[$m1]= $BA[$m1][1]*$pb[1]+$BA[$m1][2]*$pb[2]+$BA[$m1][3]*$pb[3];
        $pcnew[$m1]= $BA[$m1][1]*$pc[1]+$BA[$m1][2]*$pc[2]+$BA[$m1][3]*$pc[3];
        $pdnew[$m1]= $BA[$m1][1]*$pd[1]+$BA[$m1][2]*$pd[2]+$BA[$m1][3]*$pd[3];
        $pdnewr[$m1]= $CBA[$m1][1]*$pd[1]+$CBA[$m1][2]*$pd[2]+$CBA[$m1][3]*$pd[3];
     }

# calculando distancias entre pontos extremos a e d
  
   $dist_adnew = 0;
      $dist_adnewr= 0;
     for($m1=1; $m1<=3; $m1++){
        $da[$m1] = $pdnew[$m1] - $panew[$m1];      
        $dar[$m1]= $pdnewr[$m1] - $panew[$m1];
        $dist_adnew = $dist_adnew + $da[$m1]*$da[$m1];
        $dist_adnewr= $dist_adnewr + $dar[$m1]*$dar[$m1];
      }

    if ($dist_adnew < $dist_adnewr) { $resultado = 1;}
    else {$resultado = 0;}  
    return $resultado;
    
}