Plot()
{
 Sequence("LAT1")
 {
  Sigma_Boundary("VIa_s");
  Phase("VIa")
  {
   R_Date("65864",9018,21);
   R_Date("65863",7653,20);
   R_Date("65865",7376,19);
   R_Date("65866",5160,17);
   R_Date("65867",7740,20);
   R_Date("65868",7715,20);
   R_Date("65869",9238,22);
  };
  Sigma_Boundary("Va_s");
  Phase("Va")
  {
   R_Date("62552",6970,26);
   R_Date("54780",6998,35);
   R_Date("65862",7231,19);
   R_Date("65877",7504,20);
  };
  Sigma_Boundary("Vb_s");
  Phase("Vb")
  {
   R_Date("54765",7096,34);
   R_Date("54771",7058,33);
   R_Date("62553",7008,26);
   R_Date("65878",6903,19);
  };
  Sigma_Boundary("IV_s");
  Phase("IV")
  {
   R_Date("54775",6925,32);
   R_Date("54768",6997,32);
  };
  Sigma_Boundary("III_s");
  Phase("III")
  {
   R_Date("61861",6929,27);
   R_Date("61857",6956,27);
   R_Date("54779",5510,30);
   R_Date("61851",6918,22);
   R_Date("54774",7041,32);
   R_Date("61854",6984,26);
   R_Date("54781",6871,35);
  };
  Sigma_Boundary("III_e");
 };
 Sequence("LAT2")
 {
  Sigma_Boundary("=VIa_s");
  Phase("VIb")
  {
   R_Date("65879",6932,19);
   R_Date("65874",7068,20);
  };
  Sigma_Boundary("=Vb_s");
 };
 Sequence("LAT3")
 {
  Sigma_Boundary("=VIa_s");
  Phase("VIc")
  {
   R_Date("62551",9917,31);
  };
  Sigma_Boundary("Vc_s");
  Phase("Vc")
  {
   R_Date("62546",7078,27);
   R_Date("62547",7081,26);
   R_Date("54769",7229,33);
   R_Date("62548",7375,27);
   R_Date("62549",7112,26);
  };
  Sigma_Boundary("=IV_s");
  Phase("Vd")
  {
   R_Date("61846",6908,21);
  };
  Sigma_Boundary("=III_s");
 };
};
