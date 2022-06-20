while(<STDIN>) {
  s/agrad\((\d), (\w+)\[pp\]\)/agrad_\1_\2\[pp\]/g;
  s/grad2\((\d), (\d), (\w+)\[pp\]\)/grad2_\1_\2_\3\[pp\]/g;
  s/grad\((\d), (\w+)\[pp\]\)/grad_\1_\2\[pp\]/g;
  print $_;
}
