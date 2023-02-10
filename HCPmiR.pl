#!/usr/bin/perl -w

#Design: Chao Sun
#Program: Chao Sun
#Bug report: sunchao@caf.ac.cn

use strict;
use warnings;
use Getopt::Long;

my $usage=<<"USAGE";

Program: HCPmiR V1.0
Purpose: Identifying High-Confidence Plant miRNAs Using Multiple sRNA-Seq Data
Compiled Date: 2022-12-30

Usage: perl HCPmiR.pl -r genome -s sRNA_list -t TPM -d set_number -o out_prefix [options]*

Prerequiried tools for running: DNApi, cutadapt, Bowtie, mireap

Command-line Option
       -r    [required] Fasta format of the reference genome sequence.
       -s    [required] File list of the sRNA-Seq data.
             The suffixes of the sRNA-Seq data should be '.fastq'.
             The file list and the sRNA-Seq data should be placed in the same directory. 
       -t    [required] TPM (transcripts per million) cutoff of miRNA identification. 
       -d    [required] Set number cutoff of miRNA identification.
       -o    [required] The prefix name of the output file.
       -p    [optional] The number of threads used for Bowtie alignment.
             Default: 5.
       -m    [optional] Maximal mapping number of sRNAs on the genome.
             Default: 30.  
       -h    [Optional] Print this usage message.

Example:
       perl HCPmiR.pl -r A.genome.fa -s sRNA-Seq.list -t 2 -d 3 -o tot

USAGE

##########################Initializing input data#####################
my ($gen,$sRNA,$tpm,$set,$out,$p,$m,$help);

GetOptions(
       "r=s"=>\$gen,
       "s=s"=>\$sRNA,
       "t=i"=>\$tpm,
       "d=i"=>\$set,
       "o=s"=>\$out, 
       "p=i"=>\$p,
       "m=i"=>\$m,
       "h=s"=>\$help
);

die $usage if ($help);
die $usage if (! defined $gen);
die $usage if (! defined $sRNA);
die $usage if (! defined $tpm);
die $usage if (! defined $set);
die $usage if (! defined $out);

if(! defined $p){
   $p=5;
}
if(! defined $m){
   $m=30;
}

########Predicting and trimming adapter of sRNA-Seq data##############
my @id;
my $f=0;
my %list;
my $ada;
my $name;

open S,"$sRNA";

while(<S>){
    chomp;
    if($_=~/(\S+)\.fastq/){
       push @id,$1;
       $f ++;
       $list{$f}=$1;
    }
}
close S;

for $name(@id){
    system ("dnapi.py $name.fastq >$name.3adapt");
    open A,"$name.3adapt";

    while(<A>){
        chomp;
        $ada=$_;
    }
    close A;

    system ("cutadapt -a $ada $name.fastq -o $name.trimed.fastq");
}

##############Extracting sRNA-Seq data for mapping####################
my $file;
my $flag=1;
my $cou=1;
my %hash;
my $word;
my %s_seq;
my $num;
my $t;
my $sid;
my %read_num;

for $name(@id){
    $file = "$name.trimed.fastq";
    open IN,"$file";
    open OUT,">$name.filtered.fa";

    while(<IN>){
        chomp;
        if($flag % 4==2){
           if(length $_ >= 18 && length $_ <= 28){
              $hash{$_} ++;
           }
        }
        $flag ++;
    }

    for $word(sort keys %hash){
        $s_seq{$hash{$word}}{$word}='';
    }

    for $num(sort {$b<=>$a} keys %s_seq){
        for $t(sort {$a cmp $b} keys %{$s_seq{$num}}){
            $sid=0 x (10-length$cou).$cou;
            print OUT ">$name\_tag$sid\t$num\n$t\n";
            $read_num{$name} +=$num;
            $cou ++;
        }
    }

    close IN;
    close OUT;
    
    %hash=();
    %s_seq=();
    $flag=1;
    $cou=1;
}

#################Genome mapping and mireap############################
my $f_file;
my @cuts;
my %map;
my $start;
my $last;
my %a_pos;
my $tag;
my $m_pos;

system ("bowtie-build -f $gen GEN");

for $name(@id){
    $f_file ="$name.filtered.fa";
    system ("bowtie -a -v 0 -p $p GEN -f $f_file >$name.aln");
  
    open ALN,"$name.aln";
    open ALNOUT,">$name.aln.txt";

    while(<ALN>){
        chomp;
        @cuts=split(/\t/);
        $map{$cuts[0]} ++;
        $start=$cuts[4]+1;
        $last=$cuts[4]+(length$cuts[5]);
        $m_pos=join("\t",@cuts[0,3],$start,$last,$cuts[2]);
        push @{$a_pos{$cuts[0]}},$m_pos;
    }

    for $tag(sort keys %map){
        if($map{$tag}<=$m){
            for (@{$a_pos{$tag}}){
                print ALNOUT "$_\n";
            }
        }
    }

    close ALN;
    close ALNOUT;

    system ("mireap.pl -i $f_file -m $name.aln.txt -r $gen -o $name -t $name -u $m -d 300 -p 16");
    system ("mv $name.3adapt $name.filtered.fa $name"); 
}

######################Filitering expression values####################
my @info;
my $count;
my $pre;
my $pinfo;
my %rec;
my %MIR_pos;
my %mir;

for $name(@id){
    open R,"$name/mireap-$name.gff";
    open W,">$name/mireap-$name.TPM$tpm.filtered.gff";
    open P,">$name.TPM$tpm.pos";

    while(<R>){
        chomp;
        @info=split(/\t/);
        if($info[2] eq "precursor"){
           if($info[8]=~/ID=(\S+);(Count=)(\d+)(\S+)/){
              $count=$3*1000000/$read_num{$name};
              $info[8]="ID=$1;"."TPM=$count".$4;
              $pre=join("\t",@info);
              $MIR_pos{$1}=join("\t",@info[0,3,4,6]);
              $pinfo="$1\t$MIR_pos{$1}";
           }
        }  else {
           if($info[8]=~/(\S+)(Count=)(\d+)(\S+)/){
              $count=$3*1000000/$read_num{$name};
              $info[8]=$1."TPM=$count".$4;
              if($count >= $tpm){
                 if(! exists $rec{$pre}){
                      print W $pre,"\n";
                      print P $pinfo,"\n";
                      $rec{$pre}='';    
                 }
                 if($info[8]=~/ID=(\S+);P\S+;T\S+;Seq=(\S+)/){
                    $mir{$1}=join(";",$1,@info[0,3,4,6],$2,$count);
                 }
                 print W join("\t",@info),"\n";
              }
          }
       }
    }
    close R;
    close W;
    close P;
}

#############Integrating results of different data sets###############
my $tot='';
my $mer=1;

for $name(@id){
    system ("sort -k 2,2 -k 3,3n -k 4,4n $name.TPM$tpm.pos >$name.TPM$tpm.pos.sorted");
    &add_suff ($name,$tpm);
    &integ ($name,$tpm,\%MIR_pos,$mer);
    $tot .="$name.TPM$tpm.pos.sorted.add5.list ";
}
chop($tot);

$mer=2;
$name=$out;
system ("cat $tot >$name.TPM$tpm.pos");
system ("sort -k 2,2 -k 3,3n -k 4,4n $name.TPM$tpm.pos >$name.TPM$tpm.pos.sorted");

&add_suff ($name,$tpm);
&integ ($name,$tpm,\%MIR_pos,$mer);
$set .="sets";
&supp_miR ($name,$tpm,$set,$f,\%list,\%mir);
&check ($name,$tpm,$set);

system ("rm -rf *.ebwt *trimed.fastq *pos *pos.sort* *.aln *aln.txt");

################################Subroutine############################
sub add_suff {
    my ($name,$tpm)=@_;
    my $add=join("\t",1,1,1,1,1);

    open POS,"$name.TPM$tpm.pos.sorted";
    open POS_ADD,">$name.TPM$tpm.pos.sorted.add5";

    while(<POS>){
        chomp;
        print POS_ADD $_,"\n";
    }
    print POS_ADD $add,"\n";

    close POS;
    close POS_ADD;
}

sub integ {
    my ($name,$tpm,$MIR_pos,$mer)=@_;
    my $chr;
    my $beg;
    my $end;
    my $str;
    my $MIR;
    my $temp;
    my $number;
    my $str_num;
    my $M_out;
    my $rep; 
    my @arr;
    my @a;
    my @MIR_id;
    my @m;
    my @p_str;
    my %pos;
    my %MIR_cou;
    my %str_info;
 
    open O,"$name.TPM$tpm.pos.sorted.add5"; 
    open T,">$name.TPM$tpm.pos.sorted.add5.list";

    $temp=0;
    while(<O>){
        chomp;
        @arr=split;
        if($temp==0){
           $chr=$arr[1];
           $beg=$arr[2];
           $end=$arr[3];
           $str="$arr[4],";
           $MIR="$arr[0],";
           $temp=1;
           next;
        }

        if($arr[1] eq $chr){
           if($arr[2] > $end){
              chop($MIR);
              chop($str);
              push @{$pos{$chr}},"$beg\t$end\t$MIR\t$str";
              $beg=$arr[2];
              $end=$arr[3];
              $str="$arr[4],";
              $MIR="$arr[0],";
           }  else {
              if($arr[3] > $end){
                 $end=$arr[3];
                 $str .="$arr[4],";
                 $MIR .="$arr[0],";
              }  else {
                 $str .="$arr[4],";
                 $MIR .="$arr[0],";
              }
           }
        }  else {
             chop($MIR);
             chop($str); 
             push @{$pos{$chr}},join("\t",$beg,$end,$MIR,$str);
             $chr=$arr[1];
             $beg=$arr[2];
             $end=$arr[3];
             $str="$arr[4],";
             $MIR="$arr[0],";
        }
    }

    for my $word(sort keys %pos){
        for my $z (@{$pos{$word}}){
            @m=split(/\t/,$z);
            $M_out=join("\t",$m[2],$word,@m[0,1,3]);
            if($mer==1){
               $rep=$m[2]=~s/,/,/g;
               if($rep > 0){
                   @a=split(/,/,$m[2]);
                   print T "$a[0]\t",$MIR_pos->{$a[0]},"\n";
               } else {
                   print T "$M_out\n";
               }
            }  else {
               @p_str=split(/,/,$m[3]);
               for (@p_str){
                    $str_info{$_}='';
               } 
               $str_num=keys %str_info;
               if($str_num==1){ 
                   @MIR_id=split(/,/,$m[2]);
                   $M_out=join("\t",$m[2],$word,@m[0,1],$p_str[0]);
                   for(@MIR_id){
                       if($_=~/^(\S+)-m\d+/){
                          $MIR_cou{$1} ++;
                       }  
                   }           
                   $number=keys %{MIR_cou};
                   if($number >= $set){
                      print T $M_out,"\n";
                   }
               }
               %str_info=();
               %MIR_cou=();
            }
        }
    }
    close O;
    close T;
}

sub supp_miR {
    my ($name,$tpm,$set,$f,$list,$mir)=@_; 
    my @com;
    my @mi;
    my $o='';
    my $mid;
    my %tmp;
    my $mirna;

    open Q,"$name.TPM$tpm.pos.sorted.add5.list";
    open Y,">$name.TPM$tpm.$set.miRNA.pos";

    while(<Q>){
        chomp;
        @com=split(/\t/);
        print Y join("\t",@com),"\t";
        @mi=split(/,/,$com[0]);
        for $mid(@mi){
            if($mid=~/^(\S+)(-m\S+)/){
               $tmp{$1}=$2;
            }
        }
        for (my $j=1;$j<=$f;$j++){ 
            if(exists $tmp{$list->{$j}}){
               $mirna=$list->{$j}.$tmp{$list->{$j}};
               if(exists $mir{"$mirna-5p"}){
                  $o .=$mir{"$mirna-5p"}."\t";
               }  else {
                  $o .="-\t";
               }
               if(exists $mir{"$mirna-3p"}){
                  $o .=$mir{"$mirna-3p"}."\t";
               }  else {
                  $o .="-\t";
               }
             }  else {
                $o .="-\t-\t";
             }
         }
         chop($o);
         print Y "$o\n";
         %tmp=();
         $o='';
    }
    close Q;
    close Y;
}

sub check {
    my ($name,$tpm,$set)=@_;
    my @final;
    my @fi;
    my $f;
    my %five;
    my @left;
    my @th;
    my $h;
    my @mir_five;
    my $le_out;
    my @mir_three;
    my $ri_out;
    my %three;
    my @right;

    open C,"$name.TPM$tpm.$set.miRNA.pos";
    open H,">$name.TPM$tpm.$set.miRNA.list";

    while(<C>){
        chomp;
        @final=split(/\t/);
        print H join("\t",@final[0..4]),"\t";
        for(my $k=5;$k<=$#final;$k+=2){
            if($final[$k] ne "-"){ 
               @fi=split(/;/,$final[$k]);            
               $f=join(";",@fi[2..5]).";"; 
               $five{$fi[6]}{$f}='';
               push @left,$fi[6];
            }  else {
               push @left,"-";
            }
        }
        for(my $l=6;$l<=$#final;$l+=2){
            if($final[$l] ne "-"){ 
               @th=split(/;/,$final[$l]);            
               $h=join(";",@th[2..5]).";"; 
               $three{$th[6]}{$h}='';
               push @right,$th[6];
            }  else {
               push @right,"-";
            }
        }
        if(%five){
           @mir_five=sort {$b<=>$a} keys %five;
           for my $le(sort keys %{$five{$mir_five[0]}}){
               $le_out="$le\_".join(";",@left);
           }
        }  else {
           $le_out="NA";
        }
        if(%three){
           @mir_three=sort {$b<=>$a} keys %three;
           for my $ri(sort keys %{$three{$mir_three[0]}}){
               $ri_out="$ri\_".join(";",@right);
           }
        }  else  {
           $ri_out="NA";
        }

        print H "$le_out\t$ri_out\n";

        %five=();
        %three=();
        @left=();
        @right=();
        $le_out='';
        $ri_out='';
    }
    close C;
    close H; 
}
