\ Model amplgurobi
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x571d699d102a4c21
Minimize
  0 Assign['SML']____ + 0 Assign['JAE']____ + 0 Assign['FTR']____
   + 0 Assign['JLN']____ + 0 Assign['SDN']____ + 0 Assign['BLE']____
   + 0 Assign['ACI']____ + 0 Assign['CWU']____ + 0 Assign['JWR']____
   + 0 Assign['MBT']____ + 0 Assign['SCR']____ + 0 Assign['RCR']____
   + 0 Assign['RVY']____ + 0 Assign['PST']____ + 0 Assign['VMA']____
   + 0 Assign['RHA']____ + 0 Assign['FEE']____ + 0 Assign['RGE']____
   + 0 Assign['PSA']____ + 0 Assign['JMA']____ + 0 Assign['PSO']____
   + 0 Assign['PSN']____ + 0 Assign['JEN']____ - MinInGrp____
   + MaxInGrp____ - MinType['dept','NNE']____ - MinType['dept','NNW']____
   - MinType['dept','WSW']____ - MinType['dept','SSE']____
   - MinType['loc','Springfield']____ - MinType['loc','Urbana']____
   - MinType['loc','Peoria']____ - MinType['loc','Carbondale']____
   - MinType['loc','Joliet']____ - MinType['rate','A']____
   - MinType['rate','B']____ - MinType['title','Assistant']____
   - MinType['title','Adjunct']____ - MinType['title','Deputy']____
   + MaxType['dept','NNE']____ + MaxType['dept','NNW']____
   + MaxType['dept','WSW']____ + MaxType['dept','SSE']____
   + MaxType['loc','Springfield']____ + MaxType['loc','Urbana']____
   + MaxType['loc','Peoria']____ + MaxType['loc','Carbondale']____
   + MaxType['loc','Joliet']____ + MaxType['rate','A']____
   + MaxType['rate','B']____ + MaxType['title','Assistant']____
   + MaxType['title','Adjunct']____ + MaxType['title','Deputy']____
   + 0 MinTypeDefn[1,'dept','WSW']__3_____
   + 0 MinTypeDefn[1,'dept','WSW']__4_____
   + 0 MinTypeDefn[1,'dept','WSW']__5_____
   + 0 MinTypeDefn[1,'loc','Urbana']__5_____
   + 0 MinTypeDefn[1,'loc','Joliet']__5_____
   + 0 MinTypeDefn[2,'dept','WSW']__5_____
   + 0 MinTypeDefn[2,'loc','Urbana']__5_____
   + 0 MinTypeDefn[2,'loc','Joliet']__5_____
   + 0 MinTypeDefn[3,'dept','WSW']__5_____
   + 0 MinTypeDefn[3,'loc','Urbana']__5_____
   + 0 MinTypeDefn[3,'loc','Joliet']__5_____
   + 0 MinTypeDefn[4,'dept','WSW']__5_____
   + 0 MinTypeDefn[4,'loc','Urbana']__5_____
   + 0 MinTypeDefn[4,'loc','Joliet']__5_____
   + 0 MinTypeDefn[5,'dept','WSW']__5_____
   + 0 MinTypeDefn[5,'loc','Urbana']__5_____
   + 0 MinTypeDefn[5,'loc','Joliet']__5_____
Subject To
 MinInGrpDefn[1]__3____: MinInGrp____ - MinInGrpDefn[1]______ <= 0
 MinInGrpDefn[2]__3____: MinInGrp____ - MinInGrpDefn[2]______ <= 0
 MinInGrpDefn[3]__3____: MinInGrp____ - MinInGrpDefn[3]______ <= 0
 MinInGrpDefn[4]__3____: MinInGrp____ - MinInGrpDefn[4]______ <= 0
 MinInGrpDefn[5]__3____: MinInGrp____ - MinInGrpDefn[5]______ <= 0
 MinTypeDefn[1,'dept','NNE']__3____: MinType['dept','NNE']____
   - MinTypeDefn[1,'dept','NNE']______ <= 0
 MinTypeDefn[1,'dept','NNW']__3____: MinType['dept','NNW']____
   - MinTypeDefn[1,'dept','NNW']______ <= 0
 MinTypeDefn[1,'dept','SSE']__3____: MinType['dept','SSE']____
   - MinTypeDefn[1,'dept','SSE']______ <= 0
 MinTypeDefn[1,'loc','Springfield']__3____:
   MinType['loc','Springfield']____
   - MinTypeDefn[1,'loc','Springfield']______ <= 0
 MinTypeDefn[1,'loc','Peoria']__3____: MinType['loc','Peoria']____
   - MinTypeDefn[1,'loc','Peoria']______ <= 0
 MinTypeDefn[1,'loc','Carbondale']__3____: MinType['loc','Carbondale']____
   - MinTypeDefn[1,'loc','Carbondale']______ <= 0
 MinTypeDefn[1,'rate','A']__3____: MinType['rate','A']____
   - MinTypeDefn[1,'rate','A']______ <= 0
 MinTypeDefn[1,'rate','B']__3____: MinType['rate','B']____
   - MinTypeDefn[1,'rate','B']______ <= 0
 MinTypeDefn[1,'title','Assistant']__3____:
   MinType['title','Assistant']____
   - MinTypeDefn[1,'title','Assistant']______ <= 0
 MinTypeDefn[1,'title','Adjunct']__3____: MinType['title','Adjunct']____
   - MinTypeDefn[1,'title','Adjunct']______ <= 0
 MinTypeDefn[1,'title','Deputy']__3____: MinType['title','Deputy']____
   - MinTypeDefn[1,'title','Deputy']______ <= 0
 MinTypeDefn[2,'dept','NNE']__3____: MinType['dept','NNE']____
   - MinTypeDefn[2,'dept','NNE']______ <= 0
 MinTypeDefn[2,'dept','NNW']__3____: MinType['dept','NNW']____
   - MinTypeDefn[2,'dept','NNW']______ <= 0
 MinTypeDefn[2,'dept','SSE']__3____: MinType['dept','SSE']____
   - MinTypeDefn[2,'dept','SSE']______ <= 0
 MinTypeDefn[2,'loc','Springfield']__3____:
   MinType['loc','Springfield']____
   - MinTypeDefn[2,'loc','Springfield']______ <= 0
 MinTypeDefn[2,'loc','Peoria']__3____: MinType['loc','Peoria']____
   - MinTypeDefn[2,'loc','Peoria']______ <= 0
 MinTypeDefn[2,'loc','Carbondale']__3____: MinType['loc','Carbondale']____
   - MinTypeDefn[2,'loc','Carbondale']______ <= 0
 MinTypeDefn[2,'rate','A']__3____: MinType['rate','A']____
   - MinTypeDefn[2,'rate','A']______ <= 0
 MinTypeDefn[2,'rate','B']__3____: MinType['rate','B']____
   - MinTypeDefn[2,'rate','B']______ <= 0
 MinTypeDefn[2,'title','Assistant']__3____:
   MinType['title','Assistant']____
   - MinTypeDefn[2,'title','Assistant']______ <= 0
 MinTypeDefn[2,'title','Adjunct']__3____: MinType['title','Adjunct']____
   - MinTypeDefn[2,'title','Adjunct']______ <= 0
 MinTypeDefn[2,'title','Deputy']__3____: MinType['title','Deputy']____
   - MinTypeDefn[2,'title','Deputy']______ <= 0
 MinTypeDefn[3,'dept','NNE']__3____: MinType['dept','NNE']____
   - MinTypeDefn[3,'dept','NNE']______ <= 0
 MinTypeDefn[3,'dept','NNW']__3____: MinType['dept','NNW']____
   - MinTypeDefn[3,'dept','NNW']______ <= 0
 MinTypeDefn[3,'dept','SSE']__3____: MinType['dept','SSE']____
   - MinTypeDefn[3,'dept','SSE']______ <= 0
 MinTypeDefn[3,'loc','Springfield']__3____:
   MinType['loc','Springfield']____
   - MinTypeDefn[3,'loc','Springfield']______ <= 0
 MinTypeDefn[3,'loc','Peoria']__3____: MinType['loc','Peoria']____
   - MinTypeDefn[3,'loc','Peoria']______ <= 0
 MinTypeDefn[3,'loc','Carbondale']__3____: MinType['loc','Carbondale']____
   - MinTypeDefn[3,'loc','Carbondale']______ <= 0
 MinTypeDefn[3,'rate','A']__3____: MinType['rate','A']____
   - MinTypeDefn[3,'rate','A']______ <= 0
 MinTypeDefn[3,'rate','B']__3____: MinType['rate','B']____
   - MinTypeDefn[3,'rate','B']______ <= 0
 MinTypeDefn[3,'title','Assistant']__3____:
   MinType['title','Assistant']____
   - MinTypeDefn[3,'title','Assistant']______ <= 0
 MinTypeDefn[3,'title','Adjunct']__3____: MinType['title','Adjunct']____
   - MinTypeDefn[3,'title','Adjunct']______ <= 0
 MinTypeDefn[3,'title','Deputy']__3____: MinType['title','Deputy']____
   - MinTypeDefn[3,'title','Deputy']______ <= 0
 MinTypeDefn[4,'dept','NNE']__3____: MinType['dept','NNE']____
   - MinTypeDefn[4,'dept','NNE']______ <= 0
 MinTypeDefn[4,'dept','NNW']__3____: MinType['dept','NNW']____
   - MinTypeDefn[4,'dept','NNW']______ <= 0
 MinTypeDefn[4,'dept','SSE']__3____: MinType['dept','SSE']____
   - MinTypeDefn[4,'dept','SSE']______ <= 0
 MinTypeDefn[4,'loc','Springfield']__3____:
   MinType['loc','Springfield']____
   - MinTypeDefn[4,'loc','Springfield']______ <= 0
 MinTypeDefn[4,'loc','Peoria']__3____: MinType['loc','Peoria']____
   - MinTypeDefn[4,'loc','Peoria']______ <= 0
 MinTypeDefn[4,'loc','Carbondale']__3____: MinType['loc','Carbondale']____
   - MinTypeDefn[4,'loc','Carbondale']______ <= 0
 MinTypeDefn[4,'rate','A']__3____: MinType['rate','A']____
   - MinTypeDefn[4,'rate','A']______ <= 0
 MinTypeDefn[4,'rate','B']__3____: MinType['rate','B']____
   - MinTypeDefn[4,'rate','B']______ <= 0
 MinTypeDefn[4,'title','Assistant']__3____:
   MinType['title','Assistant']____
   - MinTypeDefn[4,'title','Assistant']______ <= 0
 MinTypeDefn[4,'title','Adjunct']__3____: MinType['title','Adjunct']____
   - MinTypeDefn[4,'title','Adjunct']______ <= 0
 MinTypeDefn[4,'title','Deputy']__3____: MinType['title','Deputy']____
   - MinTypeDefn[4,'title','Deputy']______ <= 0
 MinTypeDefn[5,'dept','NNE']__3____: MinType['dept','NNE']____
   - MinTypeDefn[5,'dept','NNE']______ <= 0
 MinTypeDefn[5,'dept','NNW']__3____: MinType['dept','NNW']____
   - MinTypeDefn[5,'dept','NNW']______ <= 0
 MinTypeDefn[5,'dept','SSE']__3____: MinType['dept','SSE']____
   - MinTypeDefn[5,'dept','SSE']______ <= 0
 MinTypeDefn[5,'loc','Springfield']__3____:
   MinType['loc','Springfield']____
   - MinTypeDefn[5,'loc','Springfield']______ <= 0
 MinTypeDefn[5,'loc','Peoria']__3____: MinType['loc','Peoria']____
   - MinTypeDefn[5,'loc','Peoria']______ <= 0
 MinTypeDefn[5,'loc','Carbondale']__3____: MinType['loc','Carbondale']____
   - MinTypeDefn[5,'loc','Carbondale']______ <= 0
 MinTypeDefn[5,'rate','A']__3____: MinType['rate','A']____
   - MinTypeDefn[5,'rate','A']______ <= 0
 MinTypeDefn[5,'rate','B']__3____: MinType['rate','B']____
   - MinTypeDefn[5,'rate','B']______ <= 0
 MinTypeDefn[5,'title','Assistant']__3____:
   MinType['title','Assistant']____
   - MinTypeDefn[5,'title','Assistant']______ <= 0
 MinTypeDefn[5,'title','Adjunct']__3____: MinType['title','Adjunct']____
   - MinTypeDefn[5,'title','Adjunct']______ <= 0
 MinTypeDefn[5,'title','Deputy']__3____: MinType['title','Deputy']____
   - MinTypeDefn[5,'title','Deputy']______ <= 0
 MinTypeDefn[1,'dept','WSW']__7_____: MinType['dept','WSW']____
   - MinTypeDefn[1,'dept','WSW']______ <= 0
 MinTypeDefn[1,'loc','Urbana']__7_____: MinType['loc','Urbana']____
   - MinTypeDefn[1,'loc','Urbana']______ <= 0
 MinTypeDefn[1,'loc','Joliet']__7_____: MinType['loc','Joliet']____
   - MinTypeDefn[1,'loc','Joliet']______ <= 0
 MinTypeDefn[2,'dept','WSW']__7_____: MinType['dept','WSW']____
   - MinTypeDefn[2,'dept','WSW']______ <= 0
 MinTypeDefn[2,'loc','Urbana']__7_____: MinType['loc','Urbana']____
   - MinTypeDefn[2,'loc','Urbana']______ <= 0
 MinTypeDefn[2,'loc','Joliet']__7_____: MinType['loc','Joliet']____
   - MinTypeDefn[2,'loc','Joliet']______ <= 0
 MinTypeDefn[3,'dept','WSW']__7_____: MinType['dept','WSW']____
   - MinTypeDefn[3,'dept','WSW']______ <= 0
 MinTypeDefn[3,'loc','Urbana']__7_____: MinType['loc','Urbana']____
   - MinTypeDefn[3,'loc','Urbana']______ <= 0
 MinTypeDefn[3,'loc','Joliet']__7_____: MinType['loc','Joliet']____
   - MinTypeDefn[3,'loc','Joliet']______ <= 0
 MinTypeDefn[4,'dept','WSW']__7_____: MinType['dept','WSW']____
   - MinTypeDefn[4,'dept','WSW']______ <= 0
 MinTypeDefn[4,'loc','Urbana']__7_____: MinType['loc','Urbana']____
   - MinTypeDefn[4,'loc','Urbana']______ <= 0
 MinTypeDefn[4,'loc','Joliet']__7_____: MinType['loc','Joliet']____
   - MinTypeDefn[4,'loc','Joliet']______ <= 0
 MinTypeDefn[5,'dept','WSW']__7_____: MinType['dept','WSW']____
   - MinTypeDefn[5,'dept','WSW']______ <= 0
 MinTypeDefn[5,'loc','Urbana']__7_____: MinType['loc','Urbana']____
   - MinTypeDefn[5,'loc','Urbana']______ <= 0
 MinTypeDefn[5,'loc','Joliet']__7_____: MinType['loc','Joliet']____
   - MinTypeDefn[5,'loc','Joliet']______ <= 0
 MinInGrpDefn[1]__2____53___: - MinInGrpDefn[1]______
   + MinTypeDefn[1,'dept','WSW']______
   + MinTypeDefn[1,'loc','Urbana']______
   + MinTypeDefn[1,'loc','Joliet']______ + MinInGrpDefn[1]__2________
   + MinInGrpDefn[1]__2____5_____ + MinInGrpDefn[1]__2____7_____
   + MinInGrpDefn[1]__2____9_____ + MinInGrpDefn[1]__2____11_____
   + MinInGrpDefn[1]__2____13_____ + MinInGrpDefn[1]__2____15_____
   + MinInGrpDefn[1]__2____17_____ + MinInGrpDefn[1]__2____19_____
   + MinInGrpDefn[1]__2____21_____ + MinInGrpDefn[1]__2____23_____
   + MinInGrpDefn[1]__2____25_____ + MinInGrpDefn[1]__2____27_____
   + MinInGrpDefn[1]__2____31_____ + MinInGrpDefn[1]__2____33_____
   + MinInGrpDefn[1]__2____35_____ + MinInGrpDefn[1]__2____37_____
   + MinInGrpDefn[1]__2____39_____ + MinInGrpDefn[1]__2____41_____
   + MinInGrpDefn[1]__2____43_____ + MinInGrpDefn[1]__2____45_____
   + MinInGrpDefn[1]__2____47_____ + MinInGrpDefn[1]__2____49_____ = 0
 MinInGrpDefn[2]__2____53___: - MinInGrpDefn[2]______
   + MinTypeDefn[2,'dept','WSW']______
   + MinTypeDefn[2,'loc','Urbana']______
   + MinTypeDefn[2,'loc','Joliet']______ + MinInGrpDefn[2]__2________
   + MinInGrpDefn[2]__2____5_____ + MinInGrpDefn[2]__2____7_____
   + MinInGrpDefn[2]__2____9_____ + MinInGrpDefn[2]__2____11_____
   + MinInGrpDefn[2]__2____13_____ + MinInGrpDefn[2]__2____15_____
   + MinInGrpDefn[2]__2____17_____ + MinInGrpDefn[2]__2____19_____
   + MinInGrpDefn[2]__2____21_____ + MinInGrpDefn[2]__2____23_____
   + MinInGrpDefn[2]__2____25_____ + MinInGrpDefn[2]__2____27_____
   + MinInGrpDefn[2]__2____31_____ + MinInGrpDefn[2]__2____33_____
   + MinInGrpDefn[2]__2____35_____ + MinInGrpDefn[2]__2____37_____
   + MinInGrpDefn[2]__2____39_____ + MinInGrpDefn[2]__2____41_____
   + MinInGrpDefn[2]__2____43_____ + MinInGrpDefn[2]__2____45_____
   + MinInGrpDefn[2]__2____47_____ + MinInGrpDefn[2]__2____49_____ = 0
 MinInGrpDefn[3]__2____53___: - MinInGrpDefn[3]______
   + MinTypeDefn[3,'dept','WSW']______
   + MinTypeDefn[3,'loc','Urbana']______
   + MinTypeDefn[3,'loc','Joliet']______ + MinInGrpDefn[3]__2________
   + MinInGrpDefn[3]__2____5_____ + MinInGrpDefn[3]__2____7_____
   + MinInGrpDefn[3]__2____9_____ + MinInGrpDefn[3]__2____11_____
   + MinInGrpDefn[3]__2____13_____ + MinInGrpDefn[3]__2____15_____
   + MinInGrpDefn[3]__2____17_____ + MinInGrpDefn[3]__2____19_____
   + MinInGrpDefn[3]__2____21_____ + MinInGrpDefn[3]__2____23_____
   + MinInGrpDefn[3]__2____25_____ + MinInGrpDefn[3]__2____27_____
   + MinInGrpDefn[3]__2____31_____ + MinInGrpDefn[3]__2____33_____
   + MinInGrpDefn[3]__2____35_____ + MinInGrpDefn[3]__2____37_____
   + MinInGrpDefn[3]__2____39_____ + MinInGrpDefn[3]__2____41_____
   + MinInGrpDefn[3]__2____43_____ + MinInGrpDefn[3]__2____45_____
   + MinInGrpDefn[3]__2____47_____ + MinInGrpDefn[3]__2____49_____ = 0
 MinInGrpDefn[4]__2____53___: - MinInGrpDefn[4]______
   + MinTypeDefn[4,'dept','WSW']______
   + MinTypeDefn[4,'loc','Urbana']______
   + MinTypeDefn[4,'loc','Joliet']______ + MinInGrpDefn[4]__2________
   + MinInGrpDefn[4]__2____5_____ + MinInGrpDefn[4]__2____7_____
   + MinInGrpDefn[4]__2____9_____ + MinInGrpDefn[4]__2____11_____
   + MinInGrpDefn[4]__2____13_____ + MinInGrpDefn[4]__2____15_____
   + MinInGrpDefn[4]__2____17_____ + MinInGrpDefn[4]__2____19_____
   + MinInGrpDefn[4]__2____21_____ + MinInGrpDefn[4]__2____23_____
   + MinInGrpDefn[4]__2____25_____ + MinInGrpDefn[4]__2____27_____
   + MinInGrpDefn[4]__2____31_____ + MinInGrpDefn[4]__2____33_____
   + MinInGrpDefn[4]__2____35_____ + MinInGrpDefn[4]__2____37_____
   + MinInGrpDefn[4]__2____39_____ + MinInGrpDefn[4]__2____41_____
   + MinInGrpDefn[4]__2____43_____ + MinInGrpDefn[4]__2____45_____
   + MinInGrpDefn[4]__2____47_____ + MinInGrpDefn[4]__2____49_____ = 0
 MinInGrpDefn[5]__2____53___: - MinInGrpDefn[5]______
   + MinTypeDefn[5,'dept','WSW']______
   + MinTypeDefn[5,'loc','Urbana']______
   + MinTypeDefn[5,'loc','Joliet']______ + MinInGrpDefn[5]__2________
   + MinInGrpDefn[5]__2____5_____ + MinInGrpDefn[5]__2____7_____
   + MinInGrpDefn[5]__2____9_____ + MinInGrpDefn[5]__2____11_____
   + MinInGrpDefn[5]__2____13_____ + MinInGrpDefn[5]__2____15_____
   + MinInGrpDefn[5]__2____17_____ + MinInGrpDefn[5]__2____19_____
   + MinInGrpDefn[5]__2____21_____ + MinInGrpDefn[5]__2____23_____
   + MinInGrpDefn[5]__2____25_____ + MinInGrpDefn[5]__2____27_____
   + MinInGrpDefn[5]__2____31_____ + MinInGrpDefn[5]__2____33_____
   + MinInGrpDefn[5]__2____35_____ + MinInGrpDefn[5]__2____37_____
   + MinInGrpDefn[5]__2____39_____ + MinInGrpDefn[5]__2____41_____
   + MinInGrpDefn[5]__2____43_____ + MinInGrpDefn[5]__2____45_____
   + MinInGrpDefn[5]__2____47_____ + MinInGrpDefn[5]__2____49_____ = 0
 MinTypeDefn[1,'dept','NNE']__2____41___:
   - MinTypeDefn[1,'dept','NNE']______
   + MinTypeDefn[1,'loc','Urbana']______
   + MinTypeDefn[1,'loc','Joliet']______ + MinInGrpDefn[1]__2________
   + MinInGrpDefn[1]__2____5_____ + MinInGrpDefn[1]__2____7_____
   + MinInGrpDefn[1]__2____9_____ + MinInGrpDefn[1]__2____11_____
   + MinInGrpDefn[1]__2____13_____ + MinInGrpDefn[1]__2____15_____
   + MinInGrpDefn[1]__2____19_____ + MinInGrpDefn[1]__2____21_____
   + MinInGrpDefn[1]__2____23_____ + MinInGrpDefn[1]__2____25_____
   + MinInGrpDefn[1]__2____27_____ + MinInGrpDefn[1]__2____31_____
   + MinInGrpDefn[1]__2____41_____ + MinInGrpDefn[1]__2____43_____
   + MinInGrpDefn[1]__2____45_____ + MinInGrpDefn[1]__2____47_____
   + MinInGrpDefn[1]__2____49_____ = 0
 MinTypeDefn[1,'dept','NNW']__2____7___:
   - MinTypeDefn[1,'dept','NNW']______ + MinInGrpDefn[1]__2____17_____
   + MinInGrpDefn[1]__2____35_____ + MinInGrpDefn[1]__2____37_____ = 0
 MinTypeDefn[1,'dept','SSE']__2____5___:
   - MinTypeDefn[1,'dept','SSE']______ + MinInGrpDefn[1]__2____33_____
   + MinInGrpDefn[1]__2____39_____ = 0
 MinTypeDefn[1,'loc','Springfield']__2____17___:
   MinTypeDefn[1,'dept','WSW']______
   - MinTypeDefn[1,'loc','Springfield']______ + MinInGrpDefn[1]__2________
   + MinInGrpDefn[1]__2____5_____ + MinInGrpDefn[1]__2____19_____
   + MinInGrpDefn[1]__2____31_____ + MinInGrpDefn[1]__2____37_____
   + MinInGrpDefn[1]__2____39_____ + MinInGrpDefn[1]__2____41_____ = 0
 MinTypeDefn[1,'loc','Peoria']__2____29___:
   - MinTypeDefn[1,'loc','Peoria']______ + MinInGrpDefn[1]__2____7_____
   + MinInGrpDefn[1]__2____9_____ + MinInGrpDefn[1]__2____11_____
   + MinInGrpDefn[1]__2____13_____ + MinInGrpDefn[1]__2____15_____
   + MinInGrpDefn[1]__2____17_____ + MinInGrpDefn[1]__2____21_____
   + MinInGrpDefn[1]__2____23_____ + MinInGrpDefn[1]__2____25_____
   + MinInGrpDefn[1]__2____27_____ + MinInGrpDefn[1]__2____33_____
   + MinInGrpDefn[1]__2____45_____ + MinInGrpDefn[1]__2____47_____
   + MinInGrpDefn[1]__2____49_____ = 0
 MinTypeDefn[1,'loc','Carbondale']__2____5___:
   - MinTypeDefn[1,'loc','Carbondale']______
   + MinInGrpDefn[1]__2____35_____ + MinInGrpDefn[1]__2____43_____ = 0
 MinTypeDefn[1,'rate','A']__2____49___: MinTypeDefn[1,'dept','WSW']______
   + MinTypeDefn[1,'loc','Urbana']______
   + MinTypeDefn[1,'loc','Joliet']______ - MinTypeDefn[1,'rate','A']______
   + MinInGrpDefn[1]__2________ + MinInGrpDefn[1]__2____5_____
   + MinInGrpDefn[1]__2____7_____ + MinInGrpDefn[1]__2____11_____
   + MinInGrpDefn[1]__2____13_____ + MinInGrpDefn[1]__2____17_____
   + MinInGrpDefn[1]__2____19_____ + MinInGrpDefn[1]__2____21_____
   + MinInGrpDefn[1]__2____23_____ + MinInGrpDefn[1]__2____25_____
   + MinInGrpDefn[1]__2____27_____ + MinInGrpDefn[1]__2____31_____
   + MinInGrpDefn[1]__2____33_____ + MinInGrpDefn[1]__2____35_____
   + MinInGrpDefn[1]__2____37_____ + MinInGrpDefn[1]__2____39_____
   + MinInGrpDefn[1]__2____41_____ + MinInGrpDefn[1]__2____43_____
   + MinInGrpDefn[1]__2____45_____ + MinInGrpDefn[1]__2____47_____
   + MinInGrpDefn[1]__2____49_____ = 0
 MinTypeDefn[1,'rate','B']__2____5___: - MinTypeDefn[1,'rate','B']______
   + MinInGrpDefn[1]__2____9_____ + MinInGrpDefn[1]__2____15_____ = 0
 MinTypeDefn[1,'title','Assistant']__2____25___:
   MinTypeDefn[1,'dept','WSW']______ + MinTypeDefn[1,'loc','Joliet']______
   - MinTypeDefn[1,'title','Assistant']______ + MinInGrpDefn[1]__2________
   + MinInGrpDefn[1]__2____5_____ + MinInGrpDefn[1]__2____9_____
   + MinInGrpDefn[1]__2____13_____ + MinInGrpDefn[1]__2____15_____
   + MinInGrpDefn[1]__2____17_____ + MinInGrpDefn[1]__2____19_____
   + MinInGrpDefn[1]__2____31_____ + MinInGrpDefn[1]__2____35_____
   + MinInGrpDefn[1]__2____41_____ = 0
 MinTypeDefn[1,'title','Adjunct']__2____23___:
   MinTypeDefn[1,'loc','Urbana']______
   - MinTypeDefn[1,'title','Adjunct']______ + MinInGrpDefn[1]__2____7_____
   + MinInGrpDefn[1]__2____21_____ + MinInGrpDefn[1]__2____23_____
   + MinInGrpDefn[1]__2____25_____ + MinInGrpDefn[1]__2____27_____
   + MinInGrpDefn[1]__2____33_____ + MinInGrpDefn[1]__2____37_____
   + MinInGrpDefn[1]__2____39_____ + MinInGrpDefn[1]__2____45_____
   + MinInGrpDefn[1]__2____47_____ = 0
 MinTypeDefn[1,'title','Deputy']__2____7___:
   - MinTypeDefn[1,'title','Deputy']______ + MinInGrpDefn[1]__2____11_____
   + MinInGrpDefn[1]__2____43_____ + MinInGrpDefn[1]__2____49_____ = 0
 MinTypeDefn[2,'dept','NNE']__2____41___:
   - MinTypeDefn[2,'dept','NNE']______
   + MinTypeDefn[2,'loc','Urbana']______
   + MinTypeDefn[2,'loc','Joliet']______ + MinInGrpDefn[2]__2________
   + MinInGrpDefn[2]__2____5_____ + MinInGrpDefn[2]__2____7_____
   + MinInGrpDefn[2]__2____9_____ + MinInGrpDefn[2]__2____11_____
   + MinInGrpDefn[2]__2____13_____ + MinInGrpDefn[2]__2____15_____
   + MinInGrpDefn[2]__2____19_____ + MinInGrpDefn[2]__2____21_____
   + MinInGrpDefn[2]__2____23_____ + MinInGrpDefn[2]__2____25_____
   + MinInGrpDefn[2]__2____27_____ + MinInGrpDefn[2]__2____31_____
   + MinInGrpDefn[2]__2____41_____ + MinInGrpDefn[2]__2____43_____
   + MinInGrpDefn[2]__2____45_____ + MinInGrpDefn[2]__2____47_____
   + MinInGrpDefn[2]__2____49_____ = 0
 MinTypeDefn[2,'dept','NNW']__2____7___:
   - MinTypeDefn[2,'dept','NNW']______ + MinInGrpDefn[2]__2____17_____
   + MinInGrpDefn[2]__2____35_____ + MinInGrpDefn[2]__2____37_____ = 0
 MinTypeDefn[2,'dept','SSE']__2____5___:
   - MinTypeDefn[2,'dept','SSE']______ + MinInGrpDefn[2]__2____33_____
   + MinInGrpDefn[2]__2____39_____ = 0
 MinTypeDefn[2,'loc','Springfield']__2____17___:
   MinTypeDefn[2,'dept','WSW']______
   - MinTypeDefn[2,'loc','Springfield']______ + MinInGrpDefn[2]__2________
   + MinInGrpDefn[2]__2____5_____ + MinInGrpDefn[2]__2____19_____
   + MinInGrpDefn[2]__2____31_____ + MinInGrpDefn[2]__2____37_____
   + MinInGrpDefn[2]__2____39_____ + MinInGrpDefn[2]__2____41_____ = 0
 MinTypeDefn[2,'loc','Peoria']__2____29___:
   - MinTypeDefn[2,'loc','Peoria']______ + MinInGrpDefn[2]__2____7_____
   + MinInGrpDefn[2]__2____9_____ + MinInGrpDefn[2]__2____11_____
   + MinInGrpDefn[2]__2____13_____ + MinInGrpDefn[2]__2____15_____
   + MinInGrpDefn[2]__2____17_____ + MinInGrpDefn[2]__2____21_____
   + MinInGrpDefn[2]__2____23_____ + MinInGrpDefn[2]__2____25_____
   + MinInGrpDefn[2]__2____27_____ + MinInGrpDefn[2]__2____33_____
   + MinInGrpDefn[2]__2____45_____ + MinInGrpDefn[2]__2____47_____
   + MinInGrpDefn[2]__2____49_____ = 0
 MinTypeDefn[2,'loc','Carbondale']__2____5___:
   - MinTypeDefn[2,'loc','Carbondale']______
   + MinInGrpDefn[2]__2____35_____ + MinInGrpDefn[2]__2____43_____ = 0
 MinTypeDefn[2,'rate','A']__2____49___: MinTypeDefn[2,'dept','WSW']______
   + MinTypeDefn[2,'loc','Urbana']______
   + MinTypeDefn[2,'loc','Joliet']______ - MinTypeDefn[2,'rate','A']______
   + MinInGrpDefn[2]__2________ + MinInGrpDefn[2]__2____5_____
   + MinInGrpDefn[2]__2____7_____ + MinInGrpDefn[2]__2____11_____
   + MinInGrpDefn[2]__2____13_____ + MinInGrpDefn[2]__2____17_____
   + MinInGrpDefn[2]__2____19_____ + MinInGrpDefn[2]__2____21_____
   + MinInGrpDefn[2]__2____23_____ + MinInGrpDefn[2]__2____25_____
   + MinInGrpDefn[2]__2____27_____ + MinInGrpDefn[2]__2____31_____
   + MinInGrpDefn[2]__2____33_____ + MinInGrpDefn[2]__2____35_____
   + MinInGrpDefn[2]__2____37_____ + MinInGrpDefn[2]__2____39_____
   + MinInGrpDefn[2]__2____41_____ + MinInGrpDefn[2]__2____43_____
   + MinInGrpDefn[2]__2____45_____ + MinInGrpDefn[2]__2____47_____
   + MinInGrpDefn[2]__2____49_____ = 0
 MinTypeDefn[2,'rate','B']__2____5___: - MinTypeDefn[2,'rate','B']______
   + MinInGrpDefn[2]__2____9_____ + MinInGrpDefn[2]__2____15_____ = 0
 MinTypeDefn[2,'title','Assistant']__2____25___:
   MinTypeDefn[2,'dept','WSW']______ + MinTypeDefn[2,'loc','Joliet']______
   - MinTypeDefn[2,'title','Assistant']______ + MinInGrpDefn[2]__2________
   + MinInGrpDefn[2]__2____5_____ + MinInGrpDefn[2]__2____9_____
   + MinInGrpDefn[2]__2____13_____ + MinInGrpDefn[2]__2____15_____
   + MinInGrpDefn[2]__2____17_____ + MinInGrpDefn[2]__2____19_____
   + MinInGrpDefn[2]__2____31_____ + MinInGrpDefn[2]__2____35_____
   + MinInGrpDefn[2]__2____41_____ = 0
 MinTypeDefn[2,'title','Adjunct']__2____23___:
   MinTypeDefn[2,'loc','Urbana']______
   - MinTypeDefn[2,'title','Adjunct']______ + MinInGrpDefn[2]__2____7_____
   + MinInGrpDefn[2]__2____21_____ + MinInGrpDefn[2]__2____23_____
   + MinInGrpDefn[2]__2____25_____ + MinInGrpDefn[2]__2____27_____
   + MinInGrpDefn[2]__2____33_____ + MinInGrpDefn[2]__2____37_____
   + MinInGrpDefn[2]__2____39_____ + MinInGrpDefn[2]__2____45_____
   + MinInGrpDefn[2]__2____47_____ = 0
 MinTypeDefn[2,'title','Deputy']__2____7___:
   - MinTypeDefn[2,'title','Deputy']______ + MinInGrpDefn[2]__2____11_____
   + MinInGrpDefn[2]__2____43_____ + MinInGrpDefn[2]__2____49_____ = 0
 MinTypeDefn[3,'dept','NNE']__2____41___:
   - MinTypeDefn[3,'dept','NNE']______
   + MinTypeDefn[3,'loc','Urbana']______
   + MinTypeDefn[3,'loc','Joliet']______ + MinInGrpDefn[3]__2________
   + MinInGrpDefn[3]__2____5_____ + MinInGrpDefn[3]__2____7_____
   + MinInGrpDefn[3]__2____9_____ + MinInGrpDefn[3]__2____11_____
   + MinInGrpDefn[3]__2____13_____ + MinInGrpDefn[3]__2____15_____
   + MinInGrpDefn[3]__2____19_____ + MinInGrpDefn[3]__2____21_____
   + MinInGrpDefn[3]__2____23_____ + MinInGrpDefn[3]__2____25_____
   + MinInGrpDefn[3]__2____27_____ + MinInGrpDefn[3]__2____31_____
   + MinInGrpDefn[3]__2____41_____ + MinInGrpDefn[3]__2____43_____
   + MinInGrpDefn[3]__2____45_____ + MinInGrpDefn[3]__2____47_____
   + MinInGrpDefn[3]__2____49_____ = 0
 MinTypeDefn[3,'dept','NNW']__2____7___:
   - MinTypeDefn[3,'dept','NNW']______ + MinInGrpDefn[3]__2____17_____
   + MinInGrpDefn[3]__2____35_____ + MinInGrpDefn[3]__2____37_____ = 0
 MinTypeDefn[3,'dept','SSE']__2____5___:
   - MinTypeDefn[3,'dept','SSE']______ + MinInGrpDefn[3]__2____33_____
   + MinInGrpDefn[3]__2____39_____ = 0
 MinTypeDefn[3,'loc','Springfield']__2____17___:
   MinTypeDefn[3,'dept','WSW']______
   - MinTypeDefn[3,'loc','Springfield']______ + MinInGrpDefn[3]__2________
   + MinInGrpDefn[3]__2____5_____ + MinInGrpDefn[3]__2____19_____
   + MinInGrpDefn[3]__2____31_____ + MinInGrpDefn[3]__2____37_____
   + MinInGrpDefn[3]__2____39_____ + MinInGrpDefn[3]__2____41_____ = 0
 MinTypeDefn[3,'loc','Peoria']__2____29___:
   - MinTypeDefn[3,'loc','Peoria']______ + MinInGrpDefn[3]__2____7_____
   + MinInGrpDefn[3]__2____9_____ + MinInGrpDefn[3]__2____11_____
   + MinInGrpDefn[3]__2____13_____ + MinInGrpDefn[3]__2____15_____
   + MinInGrpDefn[3]__2____17_____ + MinInGrpDefn[3]__2____21_____
   + MinInGrpDefn[3]__2____23_____ + MinInGrpDefn[3]__2____25_____
   + MinInGrpDefn[3]__2____27_____ + MinInGrpDefn[3]__2____33_____
   + MinInGrpDefn[3]__2____45_____ + MinInGrpDefn[3]__2____47_____
   + MinInGrpDefn[3]__2____49_____ = 0
 MinTypeDefn[3,'loc','Carbondale']__2____5___:
   - MinTypeDefn[3,'loc','Carbondale']______
   + MinInGrpDefn[3]__2____35_____ + MinInGrpDefn[3]__2____43_____ = 0
 MinTypeDefn[3,'rate','A']__2____49___: MinTypeDefn[3,'dept','WSW']______
   + MinTypeDefn[3,'loc','Urbana']______
   + MinTypeDefn[3,'loc','Joliet']______ - MinTypeDefn[3,'rate','A']______
   + MinInGrpDefn[3]__2________ + MinInGrpDefn[3]__2____5_____
   + MinInGrpDefn[3]__2____7_____ + MinInGrpDefn[3]__2____11_____
   + MinInGrpDefn[3]__2____13_____ + MinInGrpDefn[3]__2____17_____
   + MinInGrpDefn[3]__2____19_____ + MinInGrpDefn[3]__2____21_____
   + MinInGrpDefn[3]__2____23_____ + MinInGrpDefn[3]__2____25_____
   + MinInGrpDefn[3]__2____27_____ + MinInGrpDefn[3]__2____31_____
   + MinInGrpDefn[3]__2____33_____ + MinInGrpDefn[3]__2____35_____
   + MinInGrpDefn[3]__2____37_____ + MinInGrpDefn[3]__2____39_____
   + MinInGrpDefn[3]__2____41_____ + MinInGrpDefn[3]__2____43_____
   + MinInGrpDefn[3]__2____45_____ + MinInGrpDefn[3]__2____47_____
   + MinInGrpDefn[3]__2____49_____ = 0
 MinTypeDefn[3,'rate','B']__2____5___: - MinTypeDefn[3,'rate','B']______
   + MinInGrpDefn[3]__2____9_____ + MinInGrpDefn[3]__2____15_____ = 0
 MinTypeDefn[3,'title','Assistant']__2____25___:
   MinTypeDefn[3,'dept','WSW']______ + MinTypeDefn[3,'loc','Joliet']______
   - MinTypeDefn[3,'title','Assistant']______ + MinInGrpDefn[3]__2________
   + MinInGrpDefn[3]__2____5_____ + MinInGrpDefn[3]__2____9_____
   + MinInGrpDefn[3]__2____13_____ + MinInGrpDefn[3]__2____15_____
   + MinInGrpDefn[3]__2____17_____ + MinInGrpDefn[3]__2____19_____
   + MinInGrpDefn[3]__2____31_____ + MinInGrpDefn[3]__2____35_____
   + MinInGrpDefn[3]__2____41_____ = 0
 MinTypeDefn[3,'title','Adjunct']__2____23___:
   MinTypeDefn[3,'loc','Urbana']______
   - MinTypeDefn[3,'title','Adjunct']______ + MinInGrpDefn[3]__2____7_____
   + MinInGrpDefn[3]__2____21_____ + MinInGrpDefn[3]__2____23_____
   + MinInGrpDefn[3]__2____25_____ + MinInGrpDefn[3]__2____27_____
   + MinInGrpDefn[3]__2____33_____ + MinInGrpDefn[3]__2____37_____
   + MinInGrpDefn[3]__2____39_____ + MinInGrpDefn[3]__2____45_____
   + MinInGrpDefn[3]__2____47_____ = 0
 MinTypeDefn[3,'title','Deputy']__2____7___:
   - MinTypeDefn[3,'title','Deputy']______ + MinInGrpDefn[3]__2____11_____
   + MinInGrpDefn[3]__2____43_____ + MinInGrpDefn[3]__2____49_____ = 0
 MinTypeDefn[4,'dept','NNE']__2____41___:
   - MinTypeDefn[4,'dept','NNE']______
   + MinTypeDefn[4,'loc','Urbana']______
   + MinTypeDefn[4,'loc','Joliet']______ + MinInGrpDefn[4]__2________
   + MinInGrpDefn[4]__2____5_____ + MinInGrpDefn[4]__2____7_____
   + MinInGrpDefn[4]__2____9_____ + MinInGrpDefn[4]__2____11_____
   + MinInGrpDefn[4]__2____13_____ + MinInGrpDefn[4]__2____15_____
   + MinInGrpDefn[4]__2____19_____ + MinInGrpDefn[4]__2____21_____
   + MinInGrpDefn[4]__2____23_____ + MinInGrpDefn[4]__2____25_____
   + MinInGrpDefn[4]__2____27_____ + MinInGrpDefn[4]__2____31_____
   + MinInGrpDefn[4]__2____41_____ + MinInGrpDefn[4]__2____43_____
   + MinInGrpDefn[4]__2____45_____ + MinInGrpDefn[4]__2____47_____
   + MinInGrpDefn[4]__2____49_____ = 0
 MinTypeDefn[4,'dept','NNW']__2____7___:
   - MinTypeDefn[4,'dept','NNW']______ + MinInGrpDefn[4]__2____17_____
   + MinInGrpDefn[4]__2____35_____ + MinInGrpDefn[4]__2____37_____ = 0
 MinTypeDefn[4,'dept','SSE']__2____5___:
   - MinTypeDefn[4,'dept','SSE']______ + MinInGrpDefn[4]__2____33_____
   + MinInGrpDefn[4]__2____39_____ = 0
 MinTypeDefn[4,'loc','Springfield']__2____17___:
   MinTypeDefn[4,'dept','WSW']______
   - MinTypeDefn[4,'loc','Springfield']______ + MinInGrpDefn[4]__2________
   + MinInGrpDefn[4]__2____5_____ + MinInGrpDefn[4]__2____19_____
   + MinInGrpDefn[4]__2____31_____ + MinInGrpDefn[4]__2____37_____
   + MinInGrpDefn[4]__2____39_____ + MinInGrpDefn[4]__2____41_____ = 0
 MinTypeDefn[4,'loc','Peoria']__2____29___:
   - MinTypeDefn[4,'loc','Peoria']______ + MinInGrpDefn[4]__2____7_____
   + MinInGrpDefn[4]__2____9_____ + MinInGrpDefn[4]__2____11_____
   + MinInGrpDefn[4]__2____13_____ + MinInGrpDefn[4]__2____15_____
   + MinInGrpDefn[4]__2____17_____ + MinInGrpDefn[4]__2____21_____
   + MinInGrpDefn[4]__2____23_____ + MinInGrpDefn[4]__2____25_____
   + MinInGrpDefn[4]__2____27_____ + MinInGrpDefn[4]__2____33_____
   + MinInGrpDefn[4]__2____45_____ + MinInGrpDefn[4]__2____47_____
   + MinInGrpDefn[4]__2____49_____ = 0
 MinTypeDefn[4,'loc','Carbondale']__2____5___:
   - MinTypeDefn[4,'loc','Carbondale']______
   + MinInGrpDefn[4]__2____35_____ + MinInGrpDefn[4]__2____43_____ = 0
 MinTypeDefn[4,'rate','A']__2____49___: MinTypeDefn[4,'dept','WSW']______
   + MinTypeDefn[4,'loc','Urbana']______
   + MinTypeDefn[4,'loc','Joliet']______ - MinTypeDefn[4,'rate','A']______
   + MinInGrpDefn[4]__2________ + MinInGrpDefn[4]__2____5_____
   + MinInGrpDefn[4]__2____7_____ + MinInGrpDefn[4]__2____11_____
   + MinInGrpDefn[4]__2____13_____ + MinInGrpDefn[4]__2____17_____
   + MinInGrpDefn[4]__2____19_____ + MinInGrpDefn[4]__2____21_____
   + MinInGrpDefn[4]__2____23_____ + MinInGrpDefn[4]__2____25_____
   + MinInGrpDefn[4]__2____27_____ + MinInGrpDefn[4]__2____31_____
   + MinInGrpDefn[4]__2____33_____ + MinInGrpDefn[4]__2____35_____
   + MinInGrpDefn[4]__2____37_____ + MinInGrpDefn[4]__2____39_____
   + MinInGrpDefn[4]__2____41_____ + MinInGrpDefn[4]__2____43_____
   + MinInGrpDefn[4]__2____45_____ + MinInGrpDefn[4]__2____47_____
   + MinInGrpDefn[4]__2____49_____ = 0
 MinTypeDefn[4,'rate','B']__2____5___: - MinTypeDefn[4,'rate','B']______
   + MinInGrpDefn[4]__2____9_____ + MinInGrpDefn[4]__2____15_____ = 0
 MinTypeDefn[4,'title','Assistant']__2____25___:
   MinTypeDefn[4,'dept','WSW']______ + MinTypeDefn[4,'loc','Joliet']______
   - MinTypeDefn[4,'title','Assistant']______ + MinInGrpDefn[4]__2________
   + MinInGrpDefn[4]__2____5_____ + MinInGrpDefn[4]__2____9_____
   + MinInGrpDefn[4]__2____13_____ + MinInGrpDefn[4]__2____15_____
   + MinInGrpDefn[4]__2____17_____ + MinInGrpDefn[4]__2____19_____
   + MinInGrpDefn[4]__2____31_____ + MinInGrpDefn[4]__2____35_____
   + MinInGrpDefn[4]__2____41_____ = 0
 MinTypeDefn[4,'title','Adjunct']__2____23___:
   MinTypeDefn[4,'loc','Urbana']______
   - MinTypeDefn[4,'title','Adjunct']______ + MinInGrpDefn[4]__2____7_____
   + MinInGrpDefn[4]__2____21_____ + MinInGrpDefn[4]__2____23_____
   + MinInGrpDefn[4]__2____25_____ + MinInGrpDefn[4]__2____27_____
   + MinInGrpDefn[4]__2____33_____ + MinInGrpDefn[4]__2____37_____
   + MinInGrpDefn[4]__2____39_____ + MinInGrpDefn[4]__2____45_____
   + MinInGrpDefn[4]__2____47_____ = 0
 MinTypeDefn[4,'title','Deputy']__2____7___:
   - MinTypeDefn[4,'title','Deputy']______ + MinInGrpDefn[4]__2____11_____
   + MinInGrpDefn[4]__2____43_____ + MinInGrpDefn[4]__2____49_____ = 0
 MinTypeDefn[5,'dept','NNE']__2____41___:
   - MinTypeDefn[5,'dept','NNE']______
   + MinTypeDefn[5,'loc','Urbana']______
   + MinTypeDefn[5,'loc','Joliet']______ + MinInGrpDefn[5]__2________
   + MinInGrpDefn[5]__2____5_____ + MinInGrpDefn[5]__2____7_____
   + MinInGrpDefn[5]__2____9_____ + MinInGrpDefn[5]__2____11_____
   + MinInGrpDefn[5]__2____13_____ + MinInGrpDefn[5]__2____15_____
   + MinInGrpDefn[5]__2____19_____ + MinInGrpDefn[5]__2____21_____
   + MinInGrpDefn[5]__2____23_____ + MinInGrpDefn[5]__2____25_____
   + MinInGrpDefn[5]__2____27_____ + MinInGrpDefn[5]__2____31_____
   + MinInGrpDefn[5]__2____41_____ + MinInGrpDefn[5]__2____43_____
   + MinInGrpDefn[5]__2____45_____ + MinInGrpDefn[5]__2____47_____
   + MinInGrpDefn[5]__2____49_____ = 0
 MinTypeDefn[5,'dept','NNW']__2____7___:
   - MinTypeDefn[5,'dept','NNW']______ + MinInGrpDefn[5]__2____17_____
   + MinInGrpDefn[5]__2____35_____ + MinInGrpDefn[5]__2____37_____ = 0
 MinTypeDefn[5,'dept','SSE']__2____5___:
   - MinTypeDefn[5,'dept','SSE']______ + MinInGrpDefn[5]__2____33_____
   + MinInGrpDefn[5]__2____39_____ = 0
 MinTypeDefn[5,'loc','Springfield']__2____17___:
   MinTypeDefn[5,'dept','WSW']______
   - MinTypeDefn[5,'loc','Springfield']______ + MinInGrpDefn[5]__2________
   + MinInGrpDefn[5]__2____5_____ + MinInGrpDefn[5]__2____19_____
   + MinInGrpDefn[5]__2____31_____ + MinInGrpDefn[5]__2____37_____
   + MinInGrpDefn[5]__2____39_____ + MinInGrpDefn[5]__2____41_____ = 0
 MinTypeDefn[5,'loc','Peoria']__2____29___:
   - MinTypeDefn[5,'loc','Peoria']______ + MinInGrpDefn[5]__2____7_____
   + MinInGrpDefn[5]__2____9_____ + MinInGrpDefn[5]__2____11_____
   + MinInGrpDefn[5]__2____13_____ + MinInGrpDefn[5]__2____15_____
   + MinInGrpDefn[5]__2____17_____ + MinInGrpDefn[5]__2____21_____
   + MinInGrpDefn[5]__2____23_____ + MinInGrpDefn[5]__2____25_____
   + MinInGrpDefn[5]__2____27_____ + MinInGrpDefn[5]__2____33_____
   + MinInGrpDefn[5]__2____45_____ + MinInGrpDefn[5]__2____47_____
   + MinInGrpDefn[5]__2____49_____ = 0
 MinTypeDefn[5,'loc','Carbondale']__2____5___:
   - MinTypeDefn[5,'loc','Carbondale']______
   + MinInGrpDefn[5]__2____35_____ + MinInGrpDefn[5]__2____43_____ = 0
 MinTypeDefn[5,'rate','A']__2____49___: MinTypeDefn[5,'dept','WSW']______
   + MinTypeDefn[5,'loc','Urbana']______
   + MinTypeDefn[5,'loc','Joliet']______ - MinTypeDefn[5,'rate','A']______
   + MinInGrpDefn[5]__2________ + MinInGrpDefn[5]__2____5_____
   + MinInGrpDefn[5]__2____7_____ + MinInGrpDefn[5]__2____11_____
   + MinInGrpDefn[5]__2____13_____ + MinInGrpDefn[5]__2____17_____
   + MinInGrpDefn[5]__2____19_____ + MinInGrpDefn[5]__2____21_____
   + MinInGrpDefn[5]__2____23_____ + MinInGrpDefn[5]__2____25_____
   + MinInGrpDefn[5]__2____27_____ + MinInGrpDefn[5]__2____31_____
   + MinInGrpDefn[5]__2____33_____ + MinInGrpDefn[5]__2____35_____
   + MinInGrpDefn[5]__2____37_____ + MinInGrpDefn[5]__2____39_____
   + MinInGrpDefn[5]__2____41_____ + MinInGrpDefn[5]__2____43_____
   + MinInGrpDefn[5]__2____45_____ + MinInGrpDefn[5]__2____47_____
   + MinInGrpDefn[5]__2____49_____ = 0
 MinTypeDefn[5,'rate','B']__2____5___: - MinTypeDefn[5,'rate','B']______
   + MinInGrpDefn[5]__2____9_____ + MinInGrpDefn[5]__2____15_____ = 0
 MinTypeDefn[5,'title','Assistant']__2____25___:
   MinTypeDefn[5,'dept','WSW']______ + MinTypeDefn[5,'loc','Joliet']______
   - MinTypeDefn[5,'title','Assistant']______ + MinInGrpDefn[5]__2________
   + MinInGrpDefn[5]__2____5_____ + MinInGrpDefn[5]__2____9_____
   + MinInGrpDefn[5]__2____13_____ + MinInGrpDefn[5]__2____15_____
   + MinInGrpDefn[5]__2____17_____ + MinInGrpDefn[5]__2____19_____
   + MinInGrpDefn[5]__2____31_____ + MinInGrpDefn[5]__2____35_____
   + MinInGrpDefn[5]__2____41_____ = 0
 MinTypeDefn[5,'title','Adjunct']__2____23___:
   MinTypeDefn[5,'loc','Urbana']______
   - MinTypeDefn[5,'title','Adjunct']______ + MinInGrpDefn[5]__2____7_____
   + MinInGrpDefn[5]__2____21_____ + MinInGrpDefn[5]__2____23_____
   + MinInGrpDefn[5]__2____25_____ + MinInGrpDefn[5]__2____27_____
   + MinInGrpDefn[5]__2____33_____ + MinInGrpDefn[5]__2____37_____
   + MinInGrpDefn[5]__2____39_____ + MinInGrpDefn[5]__2____45_____
   + MinInGrpDefn[5]__2____47_____ = 0
 MinTypeDefn[5,'title','Deputy']__2____7___:
   - MinTypeDefn[5,'title','Deputy']______ + MinInGrpDefn[5]__2____11_____
   + MinInGrpDefn[5]__2____43_____ + MinInGrpDefn[5]__2____49_____ = 0
 MinTypeDefn[4,'loc','Joliet']__2_________:
   MinTypeDefn[1,'loc','Joliet']______
   + MinTypeDefn[2,'loc','Joliet']______
   + MinTypeDefn[3,'loc','Joliet']______
   + MinTypeDefn[4,'loc','Joliet']______
   + MinTypeDefn[5,'loc','Joliet']______ = 1
 MinTypeDefn[4,'loc','Joliet']__2_______2___: - Assign['RSN']____
   + MinTypeDefn[1,'loc','Joliet']______
   + 2 MinTypeDefn[2,'loc','Joliet']______
   + 3 MinTypeDefn[3,'loc','Joliet']______
   + 4 MinTypeDefn[4,'loc','Joliet']______
   + 5 MinTypeDefn[5,'loc','Joliet']______ = 0
 MinTypeDefn[4,'dept','WSW']__2_________: MinTypeDefn[1,'dept','WSW']______
   + MinTypeDefn[2,'dept','WSW']______ + MinTypeDefn[3,'dept','WSW']______
   + MinTypeDefn[4,'dept','WSW']______ + MinTypeDefn[5,'dept','WSW']______
   = 1
 MinTypeDefn[4,'dept','WSW']__2_______2___: - Assign['KWN']____
   + MinTypeDefn[1,'dept','WSW']______
   + 2 MinTypeDefn[2,'dept','WSW']______
   + 3 MinTypeDefn[3,'dept','WSW']______
   + 4 MinTypeDefn[4,'dept','WSW']______
   + 5 MinTypeDefn[5,'dept','WSW']______ = 0
 MinTypeDefn[4,'loc','Urbana']__2_________:
   MinTypeDefn[1,'loc','Urbana']______
   + MinTypeDefn[2,'loc','Urbana']______
   + MinTypeDefn[3,'loc','Urbana']______
   + MinTypeDefn[4,'loc','Urbana']______
   + MinTypeDefn[5,'loc','Urbana']______ = 1
 MinTypeDefn[4,'loc','Urbana']__2_______2___: - Assign['BCG']____
   + MinTypeDefn[1,'loc','Urbana']______
   + 2 MinTypeDefn[2,'loc','Urbana']______
   + 3 MinTypeDefn[3,'loc','Urbana']______
   + 4 MinTypeDefn[4,'loc','Urbana']______
   + 5 MinTypeDefn[5,'loc','Urbana']______ = 0
 MaxInGrpDefn[1]__3____: MaxInGrp____ - MinInGrpDefn[1]______ >= 0
 MaxInGrpDefn[2]__3____: MaxInGrp____ - MinInGrpDefn[2]______ >= 0
 MaxInGrpDefn[3]__3____: MaxInGrp____ - MinInGrpDefn[3]______ >= 0
 MaxInGrpDefn[4]__3____: MaxInGrp____ - MinInGrpDefn[4]______ >= 0
 MaxInGrpDefn[5]__3____: MaxInGrp____ - MinInGrpDefn[5]______ >= 0
 MaxTypeDefn[1,'dept','NNE']__3____: MaxType['dept','NNE']____
   - MinTypeDefn[1,'dept','NNE']______ >= 0
 MaxTypeDefn[1,'dept','NNW']__3____: MaxType['dept','NNW']____
   - MinTypeDefn[1,'dept','NNW']______ >= 0
 MaxTypeDefn[1,'dept','SSE']__3____: MaxType['dept','SSE']____
   - MinTypeDefn[1,'dept','SSE']______ >= 0
 MaxTypeDefn[1,'loc','Springfield']__3____:
   MaxType['loc','Springfield']____
   - MinTypeDefn[1,'loc','Springfield']______ >= 0
 MaxTypeDefn[1,'loc','Peoria']__3____: MaxType['loc','Peoria']____
   - MinTypeDefn[1,'loc','Peoria']______ >= 0
 MaxTypeDefn[1,'loc','Carbondale']__3____: MaxType['loc','Carbondale']____
   - MinTypeDefn[1,'loc','Carbondale']______ >= 0
 MaxTypeDefn[1,'rate','A']__3____: MaxType['rate','A']____
   - MinTypeDefn[1,'rate','A']______ >= 0
 MaxTypeDefn[1,'rate','B']__3____: MaxType['rate','B']____
   - MinTypeDefn[1,'rate','B']______ >= 0
 MaxTypeDefn[1,'title','Assistant']__3____:
   MaxType['title','Assistant']____
   - MinTypeDefn[1,'title','Assistant']______ >= 0
 MaxTypeDefn[1,'title','Adjunct']__3____: MaxType['title','Adjunct']____
   - MinTypeDefn[1,'title','Adjunct']______ >= 0
 MaxTypeDefn[1,'title','Deputy']__3____: MaxType['title','Deputy']____
   - MinTypeDefn[1,'title','Deputy']______ >= 0
 MaxTypeDefn[2,'dept','NNE']__3____: MaxType['dept','NNE']____
   - MinTypeDefn[2,'dept','NNE']______ >= 0
 MaxTypeDefn[2,'dept','NNW']__3____: MaxType['dept','NNW']____
   - MinTypeDefn[2,'dept','NNW']______ >= 0
 MaxTypeDefn[2,'dept','SSE']__3____: MaxType['dept','SSE']____
   - MinTypeDefn[2,'dept','SSE']______ >= 0
 MaxTypeDefn[2,'loc','Springfield']__3____:
   MaxType['loc','Springfield']____
   - MinTypeDefn[2,'loc','Springfield']______ >= 0
 MaxTypeDefn[2,'loc','Peoria']__3____: MaxType['loc','Peoria']____
   - MinTypeDefn[2,'loc','Peoria']______ >= 0
 MaxTypeDefn[2,'loc','Carbondale']__3____: MaxType['loc','Carbondale']____
   - MinTypeDefn[2,'loc','Carbondale']______ >= 0
 MaxTypeDefn[2,'rate','A']__3____: MaxType['rate','A']____
   - MinTypeDefn[2,'rate','A']______ >= 0
 MaxTypeDefn[2,'rate','B']__3____: MaxType['rate','B']____
   - MinTypeDefn[2,'rate','B']______ >= 0
 MaxTypeDefn[2,'title','Assistant']__3____:
   MaxType['title','Assistant']____
   - MinTypeDefn[2,'title','Assistant']______ >= 0
 MaxTypeDefn[2,'title','Adjunct']__3____: MaxType['title','Adjunct']____
   - MinTypeDefn[2,'title','Adjunct']______ >= 0
 MaxTypeDefn[2,'title','Deputy']__3____: MaxType['title','Deputy']____
   - MinTypeDefn[2,'title','Deputy']______ >= 0
 MaxTypeDefn[3,'dept','NNE']__3____: MaxType['dept','NNE']____
   - MinTypeDefn[3,'dept','NNE']______ >= 0
 MaxTypeDefn[3,'dept','NNW']__3____: MaxType['dept','NNW']____
   - MinTypeDefn[3,'dept','NNW']______ >= 0
 MaxTypeDefn[3,'dept','SSE']__3____: MaxType['dept','SSE']____
   - MinTypeDefn[3,'dept','SSE']______ >= 0
 MaxTypeDefn[3,'loc','Springfield']__3____:
   MaxType['loc','Springfield']____
   - MinTypeDefn[3,'loc','Springfield']______ >= 0
 MaxTypeDefn[3,'loc','Peoria']__3____: MaxType['loc','Peoria']____
   - MinTypeDefn[3,'loc','Peoria']______ >= 0
 MaxTypeDefn[3,'loc','Carbondale']__3____: MaxType['loc','Carbondale']____
   - MinTypeDefn[3,'loc','Carbondale']______ >= 0
 MaxTypeDefn[3,'rate','A']__3____: MaxType['rate','A']____
   - MinTypeDefn[3,'rate','A']______ >= 0
 MaxTypeDefn[3,'rate','B']__3____: MaxType['rate','B']____
   - MinTypeDefn[3,'rate','B']______ >= 0
 MaxTypeDefn[3,'title','Assistant']__3____:
   MaxType['title','Assistant']____
   - MinTypeDefn[3,'title','Assistant']______ >= 0
 MaxTypeDefn[3,'title','Adjunct']__3____: MaxType['title','Adjunct']____
   - MinTypeDefn[3,'title','Adjunct']______ >= 0
 MaxTypeDefn[3,'title','Deputy']__3____: MaxType['title','Deputy']____
   - MinTypeDefn[3,'title','Deputy']______ >= 0
 MaxTypeDefn[4,'dept','NNE']__3____: MaxType['dept','NNE']____
   - MinTypeDefn[4,'dept','NNE']______ >= 0
 MaxTypeDefn[4,'dept','NNW']__3____: MaxType['dept','NNW']____
   - MinTypeDefn[4,'dept','NNW']______ >= 0
 MaxTypeDefn[4,'dept','SSE']__3____: MaxType['dept','SSE']____
   - MinTypeDefn[4,'dept','SSE']______ >= 0
 MaxTypeDefn[4,'loc','Springfield']__3____:
   MaxType['loc','Springfield']____
   - MinTypeDefn[4,'loc','Springfield']______ >= 0
 MaxTypeDefn[4,'loc','Peoria']__3____: MaxType['loc','Peoria']____
   - MinTypeDefn[4,'loc','Peoria']______ >= 0
 MaxTypeDefn[4,'loc','Carbondale']__3____: MaxType['loc','Carbondale']____
   - MinTypeDefn[4,'loc','Carbondale']______ >= 0
 MaxTypeDefn[4,'rate','A']__3____: MaxType['rate','A']____
   - MinTypeDefn[4,'rate','A']______ >= 0
 MaxTypeDefn[4,'rate','B']__3____: MaxType['rate','B']____
   - MinTypeDefn[4,'rate','B']______ >= 0
 MaxTypeDefn[4,'title','Assistant']__3____:
   MaxType['title','Assistant']____
   - MinTypeDefn[4,'title','Assistant']______ >= 0
 MaxTypeDefn[4,'title','Adjunct']__3____: MaxType['title','Adjunct']____
   - MinTypeDefn[4,'title','Adjunct']______ >= 0
 MaxTypeDefn[4,'title','Deputy']__3____: MaxType['title','Deputy']____
   - MinTypeDefn[4,'title','Deputy']______ >= 0
 MaxTypeDefn[5,'dept','NNE']__3____: MaxType['dept','NNE']____
   - MinTypeDefn[5,'dept','NNE']______ >= 0
 MaxTypeDefn[5,'dept','NNW']__3____: MaxType['dept','NNW']____
   - MinTypeDefn[5,'dept','NNW']______ >= 0
 MaxTypeDefn[5,'dept','SSE']__3____: MaxType['dept','SSE']____
   - MinTypeDefn[5,'dept','SSE']______ >= 0
 MaxTypeDefn[5,'loc','Springfield']__3____:
   MaxType['loc','Springfield']____
   - MinTypeDefn[5,'loc','Springfield']______ >= 0
 MaxTypeDefn[5,'loc','Peoria']__3____: MaxType['loc','Peoria']____
   - MinTypeDefn[5,'loc','Peoria']______ >= 0
 MaxTypeDefn[5,'loc','Carbondale']__3____: MaxType['loc','Carbondale']____
   - MinTypeDefn[5,'loc','Carbondale']______ >= 0
 MaxTypeDefn[5,'rate','A']__3____: MaxType['rate','A']____
   - MinTypeDefn[5,'rate','A']______ >= 0
 MaxTypeDefn[5,'rate','B']__3____: MaxType['rate','B']____
   - MinTypeDefn[5,'rate','B']______ >= 0
 MaxTypeDefn[5,'title','Assistant']__3____:
   MaxType['title','Assistant']____
   - MinTypeDefn[5,'title','Assistant']______ >= 0
 MaxTypeDefn[5,'title','Adjunct']__3____: MaxType['title','Adjunct']____
   - MinTypeDefn[5,'title','Adjunct']______ >= 0
 MaxTypeDefn[5,'title','Deputy']__3____: MaxType['title','Deputy']____
   - MinTypeDefn[5,'title','Deputy']______ >= 0
 MaxTypeDefn[1,'dept','WSW']__7_____: MaxType['dept','WSW']____
   - MinTypeDefn[1,'dept','WSW']______ >= 0
 MaxTypeDefn[1,'loc','Urbana']__7_____: MaxType['loc','Urbana']____
   - MinTypeDefn[1,'loc','Urbana']______ >= 0
 MaxTypeDefn[1,'loc','Joliet']__7_____: MaxType['loc','Joliet']____
   - MinTypeDefn[1,'loc','Joliet']______ >= 0
 MaxTypeDefn[2,'dept','WSW']__7_____: MaxType['dept','WSW']____
   - MinTypeDefn[2,'dept','WSW']______ >= 0
 MaxTypeDefn[2,'loc','Urbana']__7_____: MaxType['loc','Urbana']____
   - MinTypeDefn[2,'loc','Urbana']______ >= 0
 MaxTypeDefn[2,'loc','Joliet']__7_____: MaxType['loc','Joliet']____
   - MinTypeDefn[2,'loc','Joliet']______ >= 0
 MaxTypeDefn[3,'dept','WSW']__7_____: MaxType['dept','WSW']____
   - MinTypeDefn[3,'dept','WSW']______ >= 0
 MaxTypeDefn[3,'loc','Urbana']__7_____: MaxType['loc','Urbana']____
   - MinTypeDefn[3,'loc','Urbana']______ >= 0
 MaxTypeDefn[3,'loc','Joliet']__7_____: MaxType['loc','Joliet']____
   - MinTypeDefn[3,'loc','Joliet']______ >= 0
 MaxTypeDefn[4,'dept','WSW']__7_____: MaxType['dept','WSW']____
   - MinTypeDefn[4,'dept','WSW']______ >= 0
 MaxTypeDefn[4,'loc','Urbana']__7_____: MaxType['loc','Urbana']____
   - MinTypeDefn[4,'loc','Urbana']______ >= 0
 MaxTypeDefn[4,'loc','Joliet']__7_____: MaxType['loc','Joliet']____
   - MinTypeDefn[4,'loc','Joliet']______ >= 0
 MaxTypeDefn[5,'dept','WSW']__7_____: MaxType['dept','WSW']____
   - MinTypeDefn[5,'dept','WSW']______ >= 0
 MaxTypeDefn[5,'loc','Urbana']__7_____: MaxType['loc','Urbana']____
   - MinTypeDefn[5,'loc','Urbana']______ >= 0
 MaxTypeDefn[5,'loc','Joliet']__7_____: MaxType['loc','Joliet']____
   - MinTypeDefn[5,'loc','Joliet']______ >= 0
 MinInGrpDefn[4]__2____50____3___: MinInGrpDefn[4]__2____49_____
   + MinInGrpDefn[4]__2____50________ + MinInGrpDefn[4]__2____50____2_____
   >= 1
 MinInGrpDefn[3]__2____50____3___: MinInGrpDefn[3]__2____49_____
   + MinInGrpDefn[3]__2____50________ + MinInGrpDefn[3]__2____50____2_____
   >= 1
 MinInGrpDefn[5]__2____50____3___: MinInGrpDefn[5]__2____49_____
   + MinInGrpDefn[5]__2____50________ + MinInGrpDefn[5]__2____50____2_____
   >= 1
 MinInGrpDefn[2]__2____50____3___: MinInGrpDefn[2]__2____49_____
   + MinInGrpDefn[2]__2____50________ + MinInGrpDefn[2]__2____50____2_____
   >= 1
 MinInGrpDefn[1]__2____50____3___: MinInGrpDefn[1]__2____49_____
   + MinInGrpDefn[1]__2____50________ + MinInGrpDefn[1]__2____50____2_____
   >= 1
 MinInGrpDefn[4]__2____48____3___: MinInGrpDefn[4]__2____47_____
   + MinInGrpDefn[4]__2____48________ + MinInGrpDefn[4]__2____48____2_____
   >= 1
 MinInGrpDefn[3]__2____48____3___: MinInGrpDefn[3]__2____47_____
   + MinInGrpDefn[3]__2____48________ + MinInGrpDefn[3]__2____48____2_____
   >= 1
 MinInGrpDefn[5]__2____48____3___: MinInGrpDefn[5]__2____47_____
   + MinInGrpDefn[5]__2____48________ + MinInGrpDefn[5]__2____48____2_____
   >= 1
 MinInGrpDefn[2]__2____48____3___: MinInGrpDefn[2]__2____47_____
   + MinInGrpDefn[2]__2____48________ + MinInGrpDefn[2]__2____48____2_____
   >= 1
 MinInGrpDefn[1]__2____48____3___: MinInGrpDefn[1]__2____47_____
   + MinInGrpDefn[1]__2____48________ + MinInGrpDefn[1]__2____48____2_____
   >= 1
 MinInGrpDefn[4]__2____46____3___: MinInGrpDefn[4]__2____45_____
   + MinInGrpDefn[4]__2____46________ + MinInGrpDefn[4]__2____46____2_____
   >= 1
 MinInGrpDefn[3]__2____46____3___: MinInGrpDefn[3]__2____45_____
   + MinInGrpDefn[3]__2____46________ + MinInGrpDefn[3]__2____46____2_____
   >= 1
 MinInGrpDefn[5]__2____46____3___: MinInGrpDefn[5]__2____45_____
   + MinInGrpDefn[5]__2____46________ + MinInGrpDefn[5]__2____46____2_____
   >= 1
 MinInGrpDefn[2]__2____46____3___: MinInGrpDefn[2]__2____45_____
   + MinInGrpDefn[2]__2____46________ + MinInGrpDefn[2]__2____46____2_____
   >= 1
 MinInGrpDefn[1]__2____46____3___: MinInGrpDefn[1]__2____45_____
   + MinInGrpDefn[1]__2____46________ + MinInGrpDefn[1]__2____46____2_____
   >= 1
 MinInGrpDefn[4]__2____44____3___: MinInGrpDefn[4]__2____43_____
   + MinInGrpDefn[4]__2____44________ + MinInGrpDefn[4]__2____44____2_____
   >= 1
 MinInGrpDefn[3]__2____44____3___: MinInGrpDefn[3]__2____43_____
   + MinInGrpDefn[3]__2____44________ + MinInGrpDefn[3]__2____44____2_____
   >= 1
 MinInGrpDefn[5]__2____44____3___: MinInGrpDefn[5]__2____43_____
   + MinInGrpDefn[5]__2____44________ + MinInGrpDefn[5]__2____44____2_____
   >= 1
 MinInGrpDefn[2]__2____44____3___: MinInGrpDefn[2]__2____43_____
   + MinInGrpDefn[2]__2____44________ + MinInGrpDefn[2]__2____44____2_____
   >= 1
 MinInGrpDefn[1]__2____44____3___: MinInGrpDefn[1]__2____43_____
   + MinInGrpDefn[1]__2____44________ + MinInGrpDefn[1]__2____44____2_____
   >= 1
 MinInGrpDefn[4]__2____42____3___: MinInGrpDefn[4]__2____41_____
   + MinInGrpDefn[4]__2____42________ + MinInGrpDefn[4]__2____42____2_____
   >= 1
 MinInGrpDefn[3]__2____42____3___: MinInGrpDefn[3]__2____41_____
   + MinInGrpDefn[3]__2____42________ + MinInGrpDefn[3]__2____42____2_____
   >= 1
 MinInGrpDefn[5]__2____42____3___: MinInGrpDefn[5]__2____41_____
   + MinInGrpDefn[5]__2____42________ + MinInGrpDefn[5]__2____42____2_____
   >= 1
 MinInGrpDefn[2]__2____42____3___: MinInGrpDefn[2]__2____41_____
   + MinInGrpDefn[2]__2____42________ + MinInGrpDefn[2]__2____42____2_____
   >= 1
 MinInGrpDefn[1]__2____42____3___: MinInGrpDefn[1]__2____41_____
   + MinInGrpDefn[1]__2____42________ + MinInGrpDefn[1]__2____42____2_____
   >= 1
 MinInGrpDefn[4]__2____40____3___: MinInGrpDefn[4]__2____39_____
   + MinInGrpDefn[4]__2____40________ + MinInGrpDefn[4]__2____40____2_____
   >= 1
 MinInGrpDefn[3]__2____40____3___: MinInGrpDefn[3]__2____39_____
   + MinInGrpDefn[3]__2____40________ + MinInGrpDefn[3]__2____40____2_____
   >= 1
 MinInGrpDefn[5]__2____40____3___: MinInGrpDefn[5]__2____39_____
   + MinInGrpDefn[5]__2____40________ + MinInGrpDefn[5]__2____40____2_____
   >= 1
 MinInGrpDefn[2]__2____40____3___: MinInGrpDefn[2]__2____39_____
   + MinInGrpDefn[2]__2____40________ + MinInGrpDefn[2]__2____40____2_____
   >= 1
 MinInGrpDefn[1]__2____40____3___: MinInGrpDefn[1]__2____39_____
   + MinInGrpDefn[1]__2____40________ + MinInGrpDefn[1]__2____40____2_____
   >= 1
 MinInGrpDefn[4]__2____38____3___: MinInGrpDefn[4]__2____37_____
   + MinInGrpDefn[4]__2____38________ + MinInGrpDefn[4]__2____38____2_____
   >= 1
 MinInGrpDefn[3]__2____38____3___: MinInGrpDefn[3]__2____37_____
   + MinInGrpDefn[3]__2____38________ + MinInGrpDefn[3]__2____38____2_____
   >= 1
 MinInGrpDefn[5]__2____38____3___: MinInGrpDefn[5]__2____37_____
   + MinInGrpDefn[5]__2____38________ + MinInGrpDefn[5]__2____38____2_____
   >= 1
 MinInGrpDefn[2]__2____38____3___: MinInGrpDefn[2]__2____37_____
   + MinInGrpDefn[2]__2____38________ + MinInGrpDefn[2]__2____38____2_____
   >= 1
 MinInGrpDefn[1]__2____38____3___: MinInGrpDefn[1]__2____37_____
   + MinInGrpDefn[1]__2____38________ + MinInGrpDefn[1]__2____38____2_____
   >= 1
 MinInGrpDefn[4]__2____36____3___: MinInGrpDefn[4]__2____35_____
   + MinInGrpDefn[4]__2____36________ + MinInGrpDefn[4]__2____36____2_____
   >= 1
 MinInGrpDefn[3]__2____36____3___: MinInGrpDefn[3]__2____35_____
   + MinInGrpDefn[3]__2____36________ + MinInGrpDefn[3]__2____36____2_____
   >= 1
 MinInGrpDefn[5]__2____36____3___: MinInGrpDefn[5]__2____35_____
   + MinInGrpDefn[5]__2____36________ + MinInGrpDefn[5]__2____36____2_____
   >= 1
 MinInGrpDefn[2]__2____36____3___: MinInGrpDefn[2]__2____35_____
   + MinInGrpDefn[2]__2____36________ + MinInGrpDefn[2]__2____36____2_____
   >= 1
 MinInGrpDefn[1]__2____36____3___: MinInGrpDefn[1]__2____35_____
   + MinInGrpDefn[1]__2____36________ + MinInGrpDefn[1]__2____36____2_____
   >= 1
 MinInGrpDefn[4]__2____34____3___: MinInGrpDefn[4]__2____33_____
   + MinInGrpDefn[4]__2____34________ + MinInGrpDefn[4]__2____34____2_____
   >= 1
 MinInGrpDefn[3]__2____34____3___: MinInGrpDefn[3]__2____33_____
   + MinInGrpDefn[3]__2____34________ + MinInGrpDefn[3]__2____34____2_____
   >= 1
 MinInGrpDefn[5]__2____34____3___: MinInGrpDefn[5]__2____33_____
   + MinInGrpDefn[5]__2____34________ + MinInGrpDefn[5]__2____34____2_____
   >= 1
 MinInGrpDefn[2]__2____34____3___: MinInGrpDefn[2]__2____33_____
   + MinInGrpDefn[2]__2____34________ + MinInGrpDefn[2]__2____34____2_____
   >= 1
 MinInGrpDefn[1]__2____34____3___: MinInGrpDefn[1]__2____33_____
   + MinInGrpDefn[1]__2____34________ + MinInGrpDefn[1]__2____34____2_____
   >= 1
 MinInGrpDefn[4]__2____32____3___: MinInGrpDefn[4]__2____31_____
   + MinInGrpDefn[4]__2____32________ + MinInGrpDefn[4]__2____32____2_____
   >= 1
 MinInGrpDefn[3]__2____32____3___: MinInGrpDefn[3]__2____31_____
   + MinInGrpDefn[3]__2____32________ + MinInGrpDefn[3]__2____32____2_____
   >= 1
 MinInGrpDefn[5]__2____32____3___: MinInGrpDefn[5]__2____31_____
   + MinInGrpDefn[5]__2____32________ + MinInGrpDefn[5]__2____32____2_____
   >= 1
 MinInGrpDefn[2]__2____32____3___: MinInGrpDefn[2]__2____31_____
   + MinInGrpDefn[2]__2____32________ + MinInGrpDefn[2]__2____32____2_____
   >= 1
 MinInGrpDefn[1]__2____32____3___: MinInGrpDefn[1]__2____31_____
   + MinInGrpDefn[1]__2____32________ + MinInGrpDefn[1]__2____32____2_____
   >= 1
 MinInGrpDefn[4]__2____28____3___: MinInGrpDefn[4]__2____27_____
   + MinInGrpDefn[4]__2____28________ + MinInGrpDefn[4]__2____28____2_____
   >= 1
 MinInGrpDefn[3]__2____28____3___: MinInGrpDefn[3]__2____27_____
   + MinInGrpDefn[3]__2____28________ + MinInGrpDefn[3]__2____28____2_____
   >= 1
 MinInGrpDefn[5]__2____28____3___: MinInGrpDefn[5]__2____27_____
   + MinInGrpDefn[5]__2____28________ + MinInGrpDefn[5]__2____28____2_____
   >= 1
 MinInGrpDefn[2]__2____28____3___: MinInGrpDefn[2]__2____27_____
   + MinInGrpDefn[2]__2____28________ + MinInGrpDefn[2]__2____28____2_____
   >= 1
 MinInGrpDefn[1]__2____28____3___: MinInGrpDefn[1]__2____27_____
   + MinInGrpDefn[1]__2____28________ + MinInGrpDefn[1]__2____28____2_____
   >= 1
 MinInGrpDefn[4]__2____26____3___: MinInGrpDefn[4]__2____25_____
   + MinInGrpDefn[4]__2____26________ + MinInGrpDefn[4]__2____26____2_____
   >= 1
 MinInGrpDefn[3]__2____26____3___: MinInGrpDefn[3]__2____25_____
   + MinInGrpDefn[3]__2____26________ + MinInGrpDefn[3]__2____26____2_____
   >= 1
 MinInGrpDefn[5]__2____26____3___: MinInGrpDefn[5]__2____25_____
   + MinInGrpDefn[5]__2____26________ + MinInGrpDefn[5]__2____26____2_____
   >= 1
 MinInGrpDefn[2]__2____26____3___: MinInGrpDefn[2]__2____25_____
   + MinInGrpDefn[2]__2____26________ + MinInGrpDefn[2]__2____26____2_____
   >= 1
 MinInGrpDefn[1]__2____26____3___: MinInGrpDefn[1]__2____25_____
   + MinInGrpDefn[1]__2____26________ + MinInGrpDefn[1]__2____26____2_____
   >= 1
 MinInGrpDefn[4]__2____24____3___: MinInGrpDefn[4]__2____23_____
   + MinInGrpDefn[4]__2____24________ + MinInGrpDefn[4]__2____24____2_____
   >= 1
 MinInGrpDefn[3]__2____24____3___: MinInGrpDefn[3]__2____23_____
   + MinInGrpDefn[3]__2____24________ + MinInGrpDefn[3]__2____24____2_____
   >= 1
 MinInGrpDefn[5]__2____24____3___: MinInGrpDefn[5]__2____23_____
   + MinInGrpDefn[5]__2____24________ + MinInGrpDefn[5]__2____24____2_____
   >= 1
 MinInGrpDefn[2]__2____24____3___: MinInGrpDefn[2]__2____23_____
   + MinInGrpDefn[2]__2____24________ + MinInGrpDefn[2]__2____24____2_____
   >= 1
 MinInGrpDefn[1]__2____24____3___: MinInGrpDefn[1]__2____23_____
   + MinInGrpDefn[1]__2____24________ + MinInGrpDefn[1]__2____24____2_____
   >= 1
 MinInGrpDefn[4]__2____22____3___: MinInGrpDefn[4]__2____21_____
   + MinInGrpDefn[4]__2____22________ + MinInGrpDefn[4]__2____22____2_____
   >= 1
 MinInGrpDefn[3]__2____22____3___: MinInGrpDefn[3]__2____21_____
   + MinInGrpDefn[3]__2____22________ + MinInGrpDefn[3]__2____22____2_____
   >= 1
 MinInGrpDefn[5]__2____22____3___: MinInGrpDefn[5]__2____21_____
   + MinInGrpDefn[5]__2____22________ + MinInGrpDefn[5]__2____22____2_____
   >= 1
 MinInGrpDefn[2]__2____22____3___: MinInGrpDefn[2]__2____21_____
   + MinInGrpDefn[2]__2____22________ + MinInGrpDefn[2]__2____22____2_____
   >= 1
 MinInGrpDefn[1]__2____22____3___: MinInGrpDefn[1]__2____21_____
   + MinInGrpDefn[1]__2____22________ + MinInGrpDefn[1]__2____22____2_____
   >= 1
 MinInGrpDefn[4]__2____20____3___: MinInGrpDefn[4]__2____19_____
   + MinInGrpDefn[4]__2____20________ + MinInGrpDefn[4]__2____20____2_____
   >= 1
 MinInGrpDefn[3]__2____20____3___: MinInGrpDefn[3]__2____19_____
   + MinInGrpDefn[3]__2____20________ + MinInGrpDefn[3]__2____20____2_____
   >= 1
 MinInGrpDefn[5]__2____20____3___: MinInGrpDefn[5]__2____19_____
   + MinInGrpDefn[5]__2____20________ + MinInGrpDefn[5]__2____20____2_____
   >= 1
 MinInGrpDefn[2]__2____20____3___: MinInGrpDefn[2]__2____19_____
   + MinInGrpDefn[2]__2____20________ + MinInGrpDefn[2]__2____20____2_____
   >= 1
 MinInGrpDefn[1]__2____20____3___: MinInGrpDefn[1]__2____19_____
   + MinInGrpDefn[1]__2____20________ + MinInGrpDefn[1]__2____20____2_____
   >= 1
 MinInGrpDefn[4]__2____18____3___: MinInGrpDefn[4]__2____17_____
   + MinInGrpDefn[4]__2____18________ + MinInGrpDefn[4]__2____18____2_____
   >= 1
 MinInGrpDefn[3]__2____18____3___: MinInGrpDefn[3]__2____17_____
   + MinInGrpDefn[3]__2____18________ + MinInGrpDefn[3]__2____18____2_____
   >= 1
 MinInGrpDefn[5]__2____18____3___: MinInGrpDefn[5]__2____17_____
   + MinInGrpDefn[5]__2____18________ + MinInGrpDefn[5]__2____18____2_____
   >= 1
 MinInGrpDefn[2]__2____18____3___: MinInGrpDefn[2]__2____17_____
   + MinInGrpDefn[2]__2____18________ + MinInGrpDefn[2]__2____18____2_____
   >= 1
 MinInGrpDefn[1]__2____18____3___: MinInGrpDefn[1]__2____17_____
   + MinInGrpDefn[1]__2____18________ + MinInGrpDefn[1]__2____18____2_____
   >= 1
 MinInGrpDefn[4]__2____16____3___: MinInGrpDefn[4]__2____15_____
   + MinInGrpDefn[4]__2____16________ + MinInGrpDefn[4]__2____16____2_____
   >= 1
 MinInGrpDefn[3]__2____16____3___: MinInGrpDefn[3]__2____15_____
   + MinInGrpDefn[3]__2____16________ + MinInGrpDefn[3]__2____16____2_____
   >= 1
 MinInGrpDefn[5]__2____16____3___: MinInGrpDefn[5]__2____15_____
   + MinInGrpDefn[5]__2____16________ + MinInGrpDefn[5]__2____16____2_____
   >= 1
 MinInGrpDefn[2]__2____16____3___: MinInGrpDefn[2]__2____15_____
   + MinInGrpDefn[2]__2____16________ + MinInGrpDefn[2]__2____16____2_____
   >= 1
 MinInGrpDefn[1]__2____16____3___: MinInGrpDefn[1]__2____15_____
   + MinInGrpDefn[1]__2____16________ + MinInGrpDefn[1]__2____16____2_____
   >= 1
 MinInGrpDefn[4]__2____14____3___: MinInGrpDefn[4]__2____13_____
   + MinInGrpDefn[4]__2____14________ + MinInGrpDefn[4]__2____14____2_____
   >= 1
 MinInGrpDefn[3]__2____14____3___: MinInGrpDefn[3]__2____13_____
   + MinInGrpDefn[3]__2____14________ + MinInGrpDefn[3]__2____14____2_____
   >= 1
 MinInGrpDefn[5]__2____14____3___: MinInGrpDefn[5]__2____13_____
   + MinInGrpDefn[5]__2____14________ + MinInGrpDefn[5]__2____14____2_____
   >= 1
 MinInGrpDefn[2]__2____14____3___: MinInGrpDefn[2]__2____13_____
   + MinInGrpDefn[2]__2____14________ + MinInGrpDefn[2]__2____14____2_____
   >= 1
 MinInGrpDefn[1]__2____14____3___: MinInGrpDefn[1]__2____13_____
   + MinInGrpDefn[1]__2____14________ + MinInGrpDefn[1]__2____14____2_____
   >= 1
 MinInGrpDefn[4]__2____12____3___: MinInGrpDefn[4]__2____11_____
   + MinInGrpDefn[4]__2____12________ + MinInGrpDefn[4]__2____12____2_____
   >= 1
 MinInGrpDefn[3]__2____12____3___: MinInGrpDefn[3]__2____11_____
   + MinInGrpDefn[3]__2____12________ + MinInGrpDefn[3]__2____12____2_____
   >= 1
 MinInGrpDefn[5]__2____12____3___: MinInGrpDefn[5]__2____11_____
   + MinInGrpDefn[5]__2____12________ + MinInGrpDefn[5]__2____12____2_____
   >= 1
 MinInGrpDefn[2]__2____12____3___: MinInGrpDefn[2]__2____11_____
   + MinInGrpDefn[2]__2____12________ + MinInGrpDefn[2]__2____12____2_____
   >= 1
 MinInGrpDefn[1]__2____12____3___: MinInGrpDefn[1]__2____11_____
   + MinInGrpDefn[1]__2____12________ + MinInGrpDefn[1]__2____12____2_____
   >= 1
 MinInGrpDefn[4]__2____10____3___: MinInGrpDefn[4]__2____9_____
   + MinInGrpDefn[4]__2____10________ + MinInGrpDefn[4]__2____10____2_____
   >= 1
 MinInGrpDefn[3]__2____10____3___: MinInGrpDefn[3]__2____9_____
   + MinInGrpDefn[3]__2____10________ + MinInGrpDefn[3]__2____10____2_____
   >= 1
 MinInGrpDefn[5]__2____10____3___: MinInGrpDefn[5]__2____9_____
   + MinInGrpDefn[5]__2____10________ + MinInGrpDefn[5]__2____10____2_____
   >= 1
 MinInGrpDefn[2]__2____10____3___: MinInGrpDefn[2]__2____9_____
   + MinInGrpDefn[2]__2____10________ + MinInGrpDefn[2]__2____10____2_____
   >= 1
 MinInGrpDefn[1]__2____10____3___: MinInGrpDefn[1]__2____9_____
   + MinInGrpDefn[1]__2____10________ + MinInGrpDefn[1]__2____10____2_____
   >= 1
 MinInGrpDefn[4]__2____6____3___: MinInGrpDefn[4]__2____5_____
   + MinInGrpDefn[4]__2____6________ + MinInGrpDefn[4]__2____6____2_____
   >= 1
 MinInGrpDefn[3]__2____6____3___: MinInGrpDefn[3]__2____5_____
   + MinInGrpDefn[3]__2____6________ + MinInGrpDefn[3]__2____6____2_____
   >= 1
 MinInGrpDefn[5]__2____6____3___: MinInGrpDefn[5]__2____5_____
   + MinInGrpDefn[5]__2____6________ + MinInGrpDefn[5]__2____6____2_____
   >= 1
 MinInGrpDefn[2]__2____6____3___: MinInGrpDefn[2]__2____5_____
   + MinInGrpDefn[2]__2____6________ + MinInGrpDefn[2]__2____6____2_____
   >= 1
 MinInGrpDefn[1]__2____6____3___: MinInGrpDefn[1]__2____5_____
   + MinInGrpDefn[1]__2____6________ + MinInGrpDefn[1]__2____6____2_____
   >= 1
 MinInGrpDefn[4]__2____2____3___: MinInGrpDefn[4]__2________
   + MinInGrpDefn[4]__2____2________ + MinInGrpDefn[4]__2____2____2_____
   >= 1
 MinInGrpDefn[3]__2____2____3___: MinInGrpDefn[3]__2________
   + MinInGrpDefn[3]__2____2________ + MinInGrpDefn[3]__2____2____2_____
   >= 1
 MinInGrpDefn[5]__2____2____3___: MinInGrpDefn[5]__2________
   + MinInGrpDefn[5]__2____2________ + MinInGrpDefn[5]__2____2____2_____
   >= 1
 MinInGrpDefn[2]__2____2____3___: MinInGrpDefn[2]__2________
   + MinInGrpDefn[2]__2____2________ + MinInGrpDefn[2]__2____2____2_____
   >= 1
 MinInGrpDefn[1]__2____2____3___: MinInGrpDefn[1]__2________
   + MinInGrpDefn[1]__2____2________ + MinInGrpDefn[1]__2____2____2_____
   >= 1
 MinInGrpDefn[4]__2____8____3___: MinInGrpDefn[4]__2____7_____
   + MinInGrpDefn[4]__2____8________ + MinInGrpDefn[4]__2____8____2_____
   >= 1
 MinInGrpDefn[3]__2____8____3___: MinInGrpDefn[3]__2____7_____
   + MinInGrpDefn[3]__2____8________ + MinInGrpDefn[3]__2____8____2_____
   >= 1
 MinInGrpDefn[5]__2____8____3___: MinInGrpDefn[5]__2____7_____
   + MinInGrpDefn[5]__2____8________ + MinInGrpDefn[5]__2____8____2_____
   >= 1
 MinInGrpDefn[2]__2____8____3___: MinInGrpDefn[2]__2____7_____
   + MinInGrpDefn[2]__2____8________ + MinInGrpDefn[2]__2____8____2_____
   >= 1
 MinInGrpDefn[1]__2____8____3___: MinInGrpDefn[1]__2____7_____
   + MinInGrpDefn[1]__2____8________ + MinInGrpDefn[1]__2____8____2_____
   >= 1
 MinInGrpDefn[4]__2____50____6___: MinInGrpDefn[4]__2____49_____ = 1 ->
   Assign['JEN']____ = 4
 MinInGrpDefn[3]__2____50____6___: MinInGrpDefn[3]__2____49_____ = 1 ->
   Assign['JEN']____ = 3
 MinInGrpDefn[5]__2____50____6___: MinInGrpDefn[5]__2____49_____ = 1 ->
   Assign['JEN']____ = 5
 MinInGrpDefn[2]__2____50____6___: MinInGrpDefn[2]__2____49_____ = 1 ->
   Assign['JEN']____ = 2
 MinInGrpDefn[1]__2____50____6___: MinInGrpDefn[1]__2____49_____ = 1 ->
   Assign['JEN']____ = 1
 MinInGrpDefn[4]__2____48____6___: MinInGrpDefn[4]__2____47_____ = 1 ->
   Assign['PSN']____ = 4
 MinInGrpDefn[3]__2____48____6___: MinInGrpDefn[3]__2____47_____ = 1 ->
   Assign['PSN']____ = 3
 MinInGrpDefn[5]__2____48____6___: MinInGrpDefn[5]__2____47_____ = 1 ->
   Assign['PSN']____ = 5
 MinInGrpDefn[2]__2____48____6___: MinInGrpDefn[2]__2____47_____ = 1 ->
   Assign['PSN']____ = 2
 MinInGrpDefn[1]__2____48____6___: MinInGrpDefn[1]__2____47_____ = 1 ->
   Assign['PSN']____ = 1
 MinInGrpDefn[4]__2____46____6___: MinInGrpDefn[4]__2____45_____ = 1 ->
   Assign['PSO']____ = 4
 MinInGrpDefn[3]__2____46____6___: MinInGrpDefn[3]__2____45_____ = 1 ->
   Assign['PSO']____ = 3
 MinInGrpDefn[5]__2____46____6___: MinInGrpDefn[5]__2____45_____ = 1 ->
   Assign['PSO']____ = 5
 MinInGrpDefn[2]__2____46____6___: MinInGrpDefn[2]__2____45_____ = 1 ->
   Assign['PSO']____ = 2
 MinInGrpDefn[1]__2____46____6___: MinInGrpDefn[1]__2____45_____ = 1 ->
   Assign['PSO']____ = 1
 MinInGrpDefn[4]__2____44____6___: MinInGrpDefn[4]__2____43_____ = 1 ->
   Assign['JMA']____ = 4
 MinInGrpDefn[3]__2____44____6___: MinInGrpDefn[3]__2____43_____ = 1 ->
   Assign['JMA']____ = 3
 MinInGrpDefn[5]__2____44____6___: MinInGrpDefn[5]__2____43_____ = 1 ->
   Assign['JMA']____ = 5
 MinInGrpDefn[2]__2____44____6___: MinInGrpDefn[2]__2____43_____ = 1 ->
   Assign['JMA']____ = 2
 MinInGrpDefn[1]__2____44____6___: MinInGrpDefn[1]__2____43_____ = 1 ->
   Assign['JMA']____ = 1
 MinInGrpDefn[4]__2____42____6___: MinInGrpDefn[4]__2____41_____ = 1 ->
   Assign['PSA']____ = 4
 MinInGrpDefn[3]__2____42____6___: MinInGrpDefn[3]__2____41_____ = 1 ->
   Assign['PSA']____ = 3
 MinInGrpDefn[5]__2____42____6___: MinInGrpDefn[5]__2____41_____ = 1 ->
   Assign['PSA']____ = 5
 MinInGrpDefn[2]__2____42____6___: MinInGrpDefn[2]__2____41_____ = 1 ->
   Assign['PSA']____ = 2
 MinInGrpDefn[1]__2____42____6___: MinInGrpDefn[1]__2____41_____ = 1 ->
   Assign['PSA']____ = 1
 MinInGrpDefn[4]__2____40____6___: MinInGrpDefn[4]__2____39_____ = 1 ->
   Assign['RGE']____ = 4
 MinInGrpDefn[3]__2____40____6___: MinInGrpDefn[3]__2____39_____ = 1 ->
   Assign['RGE']____ = 3
 MinInGrpDefn[5]__2____40____6___: MinInGrpDefn[5]__2____39_____ = 1 ->
   Assign['RGE']____ = 5
 MinInGrpDefn[2]__2____40____6___: MinInGrpDefn[2]__2____39_____ = 1 ->
   Assign['RGE']____ = 2
 MinInGrpDefn[1]__2____40____6___: MinInGrpDefn[1]__2____39_____ = 1 ->
   Assign['RGE']____ = 1
 MinInGrpDefn[4]__2____38____6___: MinInGrpDefn[4]__2____37_____ = 1 ->
   Assign['FEE']____ = 4
 MinInGrpDefn[3]__2____38____6___: MinInGrpDefn[3]__2____37_____ = 1 ->
   Assign['FEE']____ = 3
 MinInGrpDefn[5]__2____38____6___: MinInGrpDefn[5]__2____37_____ = 1 ->
   Assign['FEE']____ = 5
 MinInGrpDefn[2]__2____38____6___: MinInGrpDefn[2]__2____37_____ = 1 ->
   Assign['FEE']____ = 2
 MinInGrpDefn[1]__2____38____6___: MinInGrpDefn[1]__2____37_____ = 1 ->
   Assign['FEE']____ = 1
 MinInGrpDefn[4]__2____36____6___: MinInGrpDefn[4]__2____35_____ = 1 ->
   Assign['RHA']____ = 4
 MinInGrpDefn[3]__2____36____6___: MinInGrpDefn[3]__2____35_____ = 1 ->
   Assign['RHA']____ = 3
 MinInGrpDefn[5]__2____36____6___: MinInGrpDefn[5]__2____35_____ = 1 ->
   Assign['RHA']____ = 5
 MinInGrpDefn[2]__2____36____6___: MinInGrpDefn[2]__2____35_____ = 1 ->
   Assign['RHA']____ = 2
 MinInGrpDefn[1]__2____36____6___: MinInGrpDefn[1]__2____35_____ = 1 ->
   Assign['RHA']____ = 1
 MinInGrpDefn[4]__2____34____6___: MinInGrpDefn[4]__2____33_____ = 1 ->
   Assign['VMA']____ = 4
 MinInGrpDefn[3]__2____34____6___: MinInGrpDefn[3]__2____33_____ = 1 ->
   Assign['VMA']____ = 3
 MinInGrpDefn[5]__2____34____6___: MinInGrpDefn[5]__2____33_____ = 1 ->
   Assign['VMA']____ = 5
 MinInGrpDefn[2]__2____34____6___: MinInGrpDefn[2]__2____33_____ = 1 ->
   Assign['VMA']____ = 2
 MinInGrpDefn[1]__2____34____6___: MinInGrpDefn[1]__2____33_____ = 1 ->
   Assign['VMA']____ = 1
 MinInGrpDefn[4]__2____32____6___: MinInGrpDefn[4]__2____31_____ = 1 ->
   Assign['PST']____ = 4
 MinInGrpDefn[3]__2____32____6___: MinInGrpDefn[3]__2____31_____ = 1 ->
   Assign['PST']____ = 3
 MinInGrpDefn[5]__2____32____6___: MinInGrpDefn[5]__2____31_____ = 1 ->
   Assign['PST']____ = 5
 MinInGrpDefn[2]__2____32____6___: MinInGrpDefn[2]__2____31_____ = 1 ->
   Assign['PST']____ = 2
 MinInGrpDefn[1]__2____32____6___: MinInGrpDefn[1]__2____31_____ = 1 ->
   Assign['PST']____ = 1
 MinInGrpDefn[4]__2____28____6___: MinInGrpDefn[4]__2____27_____ = 1 ->
   Assign['RVY']____ = 4
 MinInGrpDefn[3]__2____28____6___: MinInGrpDefn[3]__2____27_____ = 1 ->
   Assign['RVY']____ = 3
 MinInGrpDefn[5]__2____28____6___: MinInGrpDefn[5]__2____27_____ = 1 ->
   Assign['RVY']____ = 5
 MinInGrpDefn[2]__2____28____6___: MinInGrpDefn[2]__2____27_____ = 1 ->
   Assign['RVY']____ = 2
 MinInGrpDefn[1]__2____28____6___: MinInGrpDefn[1]__2____27_____ = 1 ->
   Assign['RVY']____ = 1
 MinInGrpDefn[4]__2____26____6___: MinInGrpDefn[4]__2____25_____ = 1 ->
   Assign['RCR']____ = 4
 MinInGrpDefn[3]__2____26____6___: MinInGrpDefn[3]__2____25_____ = 1 ->
   Assign['RCR']____ = 3
 MinInGrpDefn[5]__2____26____6___: MinInGrpDefn[5]__2____25_____ = 1 ->
   Assign['RCR']____ = 5
 MinInGrpDefn[2]__2____26____6___: MinInGrpDefn[2]__2____25_____ = 1 ->
   Assign['RCR']____ = 2
 MinInGrpDefn[1]__2____26____6___: MinInGrpDefn[1]__2____25_____ = 1 ->
   Assign['RCR']____ = 1
 MinInGrpDefn[4]__2____24____6___: MinInGrpDefn[4]__2____23_____ = 1 ->
   Assign['SCR']____ = 4
 MinInGrpDefn[3]__2____24____6___: MinInGrpDefn[3]__2____23_____ = 1 ->
   Assign['SCR']____ = 3
 MinInGrpDefn[5]__2____24____6___: MinInGrpDefn[5]__2____23_____ = 1 ->
   Assign['SCR']____ = 5
 MinInGrpDefn[2]__2____24____6___: MinInGrpDefn[2]__2____23_____ = 1 ->
   Assign['SCR']____ = 2
 MinInGrpDefn[1]__2____24____6___: MinInGrpDefn[1]__2____23_____ = 1 ->
   Assign['SCR']____ = 1
 MinInGrpDefn[4]__2____22____6___: MinInGrpDefn[4]__2____21_____ = 1 ->
   Assign['MBT']____ = 4
 MinInGrpDefn[3]__2____22____6___: MinInGrpDefn[3]__2____21_____ = 1 ->
   Assign['MBT']____ = 3
 MinInGrpDefn[5]__2____22____6___: MinInGrpDefn[5]__2____21_____ = 1 ->
   Assign['MBT']____ = 5
 MinInGrpDefn[2]__2____22____6___: MinInGrpDefn[2]__2____21_____ = 1 ->
   Assign['MBT']____ = 2
 MinInGrpDefn[1]__2____22____6___: MinInGrpDefn[1]__2____21_____ = 1 ->
   Assign['MBT']____ = 1
 MinInGrpDefn[4]__2____20____6___: MinInGrpDefn[4]__2____19_____ = 1 ->
   Assign['JWR']____ = 4
 MinInGrpDefn[3]__2____20____6___: MinInGrpDefn[3]__2____19_____ = 1 ->
   Assign['JWR']____ = 3
 MinInGrpDefn[5]__2____20____6___: MinInGrpDefn[5]__2____19_____ = 1 ->
   Assign['JWR']____ = 5
 MinInGrpDefn[2]__2____20____6___: MinInGrpDefn[2]__2____19_____ = 1 ->
   Assign['JWR']____ = 2
 MinInGrpDefn[1]__2____20____6___: MinInGrpDefn[1]__2____19_____ = 1 ->
   Assign['JWR']____ = 1
 MinInGrpDefn[4]__2____18____6___: MinInGrpDefn[4]__2____17_____ = 1 ->
   Assign['CWU']____ = 4
 MinInGrpDefn[3]__2____18____6___: MinInGrpDefn[3]__2____17_____ = 1 ->
   Assign['CWU']____ = 3
 MinInGrpDefn[5]__2____18____6___: MinInGrpDefn[5]__2____17_____ = 1 ->
   Assign['CWU']____ = 5
 MinInGrpDefn[2]__2____18____6___: MinInGrpDefn[2]__2____17_____ = 1 ->
   Assign['CWU']____ = 2
 MinInGrpDefn[1]__2____18____6___: MinInGrpDefn[1]__2____17_____ = 1 ->
   Assign['CWU']____ = 1
 MinInGrpDefn[4]__2____16____6___: MinInGrpDefn[4]__2____15_____ = 1 ->
   Assign['ACI']____ = 4
 MinInGrpDefn[3]__2____16____6___: MinInGrpDefn[3]__2____15_____ = 1 ->
   Assign['ACI']____ = 3
 MinInGrpDefn[5]__2____16____6___: MinInGrpDefn[5]__2____15_____ = 1 ->
   Assign['ACI']____ = 5
 MinInGrpDefn[2]__2____16____6___: MinInGrpDefn[2]__2____15_____ = 1 ->
   Assign['ACI']____ = 2
 MinInGrpDefn[1]__2____16____6___: MinInGrpDefn[1]__2____15_____ = 1 ->
   Assign['ACI']____ = 1
 MinInGrpDefn[4]__2____14____6___: MinInGrpDefn[4]__2____13_____ = 1 ->
   Assign['BLE']____ = 4
 MinInGrpDefn[3]__2____14____6___: MinInGrpDefn[3]__2____13_____ = 1 ->
   Assign['BLE']____ = 3
 MinInGrpDefn[5]__2____14____6___: MinInGrpDefn[5]__2____13_____ = 1 ->
   Assign['BLE']____ = 5
 MinInGrpDefn[2]__2____14____6___: MinInGrpDefn[2]__2____13_____ = 1 ->
   Assign['BLE']____ = 2
 MinInGrpDefn[1]__2____14____6___: MinInGrpDefn[1]__2____13_____ = 1 ->
   Assign['BLE']____ = 1
 MinInGrpDefn[4]__2____12____6___: MinInGrpDefn[4]__2____11_____ = 1 ->
   Assign['SDN']____ = 4
 MinInGrpDefn[3]__2____12____6___: MinInGrpDefn[3]__2____11_____ = 1 ->
   Assign['SDN']____ = 3
 MinInGrpDefn[5]__2____12____6___: MinInGrpDefn[5]__2____11_____ = 1 ->
   Assign['SDN']____ = 5
 MinInGrpDefn[2]__2____12____6___: MinInGrpDefn[2]__2____11_____ = 1 ->
   Assign['SDN']____ = 2
 MinInGrpDefn[1]__2____12____6___: MinInGrpDefn[1]__2____11_____ = 1 ->
   Assign['SDN']____ = 1
 MinInGrpDefn[4]__2____10____6___: MinInGrpDefn[4]__2____9_____ = 1 ->
   Assign['JLN']____ = 4
 MinInGrpDefn[3]__2____10____6___: MinInGrpDefn[3]__2____9_____ = 1 ->
   Assign['JLN']____ = 3
 MinInGrpDefn[5]__2____10____6___: MinInGrpDefn[5]__2____9_____ = 1 ->
   Assign['JLN']____ = 5
 MinInGrpDefn[2]__2____10____6___: MinInGrpDefn[2]__2____9_____ = 1 ->
   Assign['JLN']____ = 2
 MinInGrpDefn[1]__2____10____6___: MinInGrpDefn[1]__2____9_____ = 1 ->
   Assign['JLN']____ = 1
 MinInGrpDefn[4]__2____6____6___: MinInGrpDefn[4]__2____5_____ = 1 ->
   Assign['JAE']____ = 4
 MinInGrpDefn[3]__2____6____6___: MinInGrpDefn[3]__2____5_____ = 1 ->
   Assign['JAE']____ = 3
 MinInGrpDefn[5]__2____6____6___: MinInGrpDefn[5]__2____5_____ = 1 ->
   Assign['JAE']____ = 5
 MinInGrpDefn[2]__2____6____6___: MinInGrpDefn[2]__2____5_____ = 1 ->
   Assign['JAE']____ = 2
 MinInGrpDefn[1]__2____6____6___: MinInGrpDefn[1]__2____5_____ = 1 ->
   Assign['JAE']____ = 1
 MinInGrpDefn[4]__2____2____6___: MinInGrpDefn[4]__2________ = 1 ->
   Assign['SML']____ = 4
 MinInGrpDefn[3]__2____2____6___: MinInGrpDefn[3]__2________ = 1 ->
   Assign['SML']____ = 3
 MinInGrpDefn[5]__2____2____6___: MinInGrpDefn[5]__2________ = 1 ->
   Assign['SML']____ = 5
 MinInGrpDefn[2]__2____2____6___: MinInGrpDefn[2]__2________ = 1 ->
   Assign['SML']____ = 2
 MinInGrpDefn[1]__2____2____6___: MinInGrpDefn[1]__2________ = 1 ->
   Assign['SML']____ = 1
 MinInGrpDefn[4]__2____8____6___: MinInGrpDefn[4]__2____7_____ = 1 ->
   Assign['FTR']____ = 4
 MinInGrpDefn[3]__2____8____6___: MinInGrpDefn[3]__2____7_____ = 1 ->
   Assign['FTR']____ = 3
 MinInGrpDefn[5]__2____8____6___: MinInGrpDefn[5]__2____7_____ = 1 ->
   Assign['FTR']____ = 5
 MinInGrpDefn[2]__2____8____6___: MinInGrpDefn[2]__2____7_____ = 1 ->
   Assign['FTR']____ = 2
 MinInGrpDefn[1]__2____8____6___: MinInGrpDefn[1]__2____7_____ = 1 ->
   Assign['FTR']____ = 1
 MinInGrpDefn[4]__2____50____4___: MinInGrpDefn[4]__2____50________ = 1 ->
   Assign['JEN']____ <= 3
 MinInGrpDefn[3]__2____50____4___: MinInGrpDefn[3]__2____50________ = 1 ->
   Assign['JEN']____ <= 2
 MinInGrpDefn[5]__2____50____4___: MinInGrpDefn[5]__2____50________ = 1 ->
   Assign['JEN']____ <= 4
 MinInGrpDefn[2]__2____50____4___: MinInGrpDefn[2]__2____50________ = 1 ->
   Assign['JEN']____ <= 1
 MinInGrpDefn[1]__2____50____4___: MinInGrpDefn[1]__2____50________ = 1 ->
   Assign['JEN']____ <= 0
 MinInGrpDefn[4]__2____48____4___: MinInGrpDefn[4]__2____48________ = 1 ->
   Assign['PSN']____ <= 3
 MinInGrpDefn[3]__2____48____4___: MinInGrpDefn[3]__2____48________ = 1 ->
   Assign['PSN']____ <= 2
 MinInGrpDefn[5]__2____48____4___: MinInGrpDefn[5]__2____48________ = 1 ->
   Assign['PSN']____ <= 4
 MinInGrpDefn[2]__2____48____4___: MinInGrpDefn[2]__2____48________ = 1 ->
   Assign['PSN']____ <= 1
 MinInGrpDefn[1]__2____48____4___: MinInGrpDefn[1]__2____48________ = 1 ->
   Assign['PSN']____ <= 0
 MinInGrpDefn[4]__2____46____4___: MinInGrpDefn[4]__2____46________ = 1 ->
   Assign['PSO']____ <= 3
 MinInGrpDefn[3]__2____46____4___: MinInGrpDefn[3]__2____46________ = 1 ->
   Assign['PSO']____ <= 2
 MinInGrpDefn[5]__2____46____4___: MinInGrpDefn[5]__2____46________ = 1 ->
   Assign['PSO']____ <= 4
 MinInGrpDefn[2]__2____46____4___: MinInGrpDefn[2]__2____46________ = 1 ->
   Assign['PSO']____ <= 1
 MinInGrpDefn[1]__2____46____4___: MinInGrpDefn[1]__2____46________ = 1 ->
   Assign['PSO']____ <= 0
 MinInGrpDefn[4]__2____44____4___: MinInGrpDefn[4]__2____44________ = 1 ->
   Assign['JMA']____ <= 3
 MinInGrpDefn[3]__2____44____4___: MinInGrpDefn[3]__2____44________ = 1 ->
   Assign['JMA']____ <= 2
 MinInGrpDefn[5]__2____44____4___: MinInGrpDefn[5]__2____44________ = 1 ->
   Assign['JMA']____ <= 4
 MinInGrpDefn[2]__2____44____4___: MinInGrpDefn[2]__2____44________ = 1 ->
   Assign['JMA']____ <= 1
 MinInGrpDefn[1]__2____44____4___: MinInGrpDefn[1]__2____44________ = 1 ->
   Assign['JMA']____ <= 0
 MinInGrpDefn[4]__2____42____4___: MinInGrpDefn[4]__2____42________ = 1 ->
   Assign['PSA']____ <= 3
 MinInGrpDefn[3]__2____42____4___: MinInGrpDefn[3]__2____42________ = 1 ->
   Assign['PSA']____ <= 2
 MinInGrpDefn[5]__2____42____4___: MinInGrpDefn[5]__2____42________ = 1 ->
   Assign['PSA']____ <= 4
 MinInGrpDefn[2]__2____42____4___: MinInGrpDefn[2]__2____42________ = 1 ->
   Assign['PSA']____ <= 1
 MinInGrpDefn[1]__2____42____4___: MinInGrpDefn[1]__2____42________ = 1 ->
   Assign['PSA']____ <= 0
 MinInGrpDefn[4]__2____40____4___: MinInGrpDefn[4]__2____40________ = 1 ->
   Assign['RGE']____ <= 3
 MinInGrpDefn[3]__2____40____4___: MinInGrpDefn[3]__2____40________ = 1 ->
   Assign['RGE']____ <= 2
 MinInGrpDefn[5]__2____40____4___: MinInGrpDefn[5]__2____40________ = 1 ->
   Assign['RGE']____ <= 4
 MinInGrpDefn[2]__2____40____4___: MinInGrpDefn[2]__2____40________ = 1 ->
   Assign['RGE']____ <= 1
 MinInGrpDefn[1]__2____40____4___: MinInGrpDefn[1]__2____40________ = 1 ->
   Assign['RGE']____ <= 0
 MinInGrpDefn[4]__2____38____4___: MinInGrpDefn[4]__2____38________ = 1 ->
   Assign['FEE']____ <= 3
 MinInGrpDefn[3]__2____38____4___: MinInGrpDefn[3]__2____38________ = 1 ->
   Assign['FEE']____ <= 2
 MinInGrpDefn[5]__2____38____4___: MinInGrpDefn[5]__2____38________ = 1 ->
   Assign['FEE']____ <= 4
 MinInGrpDefn[2]__2____38____4___: MinInGrpDefn[2]__2____38________ = 1 ->
   Assign['FEE']____ <= 1
 MinInGrpDefn[1]__2____38____4___: MinInGrpDefn[1]__2____38________ = 1 ->
   Assign['FEE']____ <= 0
 MinInGrpDefn[4]__2____36____4___: MinInGrpDefn[4]__2____36________ = 1 ->
   Assign['RHA']____ <= 3
 MinInGrpDefn[3]__2____36____4___: MinInGrpDefn[3]__2____36________ = 1 ->
   Assign['RHA']____ <= 2
 MinInGrpDefn[5]__2____36____4___: MinInGrpDefn[5]__2____36________ = 1 ->
   Assign['RHA']____ <= 4
 MinInGrpDefn[2]__2____36____4___: MinInGrpDefn[2]__2____36________ = 1 ->
   Assign['RHA']____ <= 1
 MinInGrpDefn[1]__2____36____4___: MinInGrpDefn[1]__2____36________ = 1 ->
   Assign['RHA']____ <= 0
 MinInGrpDefn[4]__2____34____4___: MinInGrpDefn[4]__2____34________ = 1 ->
   Assign['VMA']____ <= 3
 MinInGrpDefn[3]__2____34____4___: MinInGrpDefn[3]__2____34________ = 1 ->
   Assign['VMA']____ <= 2
 MinInGrpDefn[5]__2____34____4___: MinInGrpDefn[5]__2____34________ = 1 ->
   Assign['VMA']____ <= 4
 MinInGrpDefn[2]__2____34____4___: MinInGrpDefn[2]__2____34________ = 1 ->
   Assign['VMA']____ <= 1
 MinInGrpDefn[1]__2____34____4___: MinInGrpDefn[1]__2____34________ = 1 ->
   Assign['VMA']____ <= 0
 MinInGrpDefn[4]__2____32____4___: MinInGrpDefn[4]__2____32________ = 1 ->
   Assign['PST']____ <= 3
 MinInGrpDefn[3]__2____32____4___: MinInGrpDefn[3]__2____32________ = 1 ->
   Assign['PST']____ <= 2
 MinInGrpDefn[5]__2____32____4___: MinInGrpDefn[5]__2____32________ = 1 ->
   Assign['PST']____ <= 4
 MinInGrpDefn[2]__2____32____4___: MinInGrpDefn[2]__2____32________ = 1 ->
   Assign['PST']____ <= 1
 MinInGrpDefn[1]__2____32____4___: MinInGrpDefn[1]__2____32________ = 1 ->
   Assign['PST']____ <= 0
 MinInGrpDefn[4]__2____28____4___: MinInGrpDefn[4]__2____28________ = 1 ->
   Assign['RVY']____ <= 3
 MinInGrpDefn[3]__2____28____4___: MinInGrpDefn[3]__2____28________ = 1 ->
   Assign['RVY']____ <= 2
 MinInGrpDefn[5]__2____28____4___: MinInGrpDefn[5]__2____28________ = 1 ->
   Assign['RVY']____ <= 4
 MinInGrpDefn[2]__2____28____4___: MinInGrpDefn[2]__2____28________ = 1 ->
   Assign['RVY']____ <= 1
 MinInGrpDefn[1]__2____28____4___: MinInGrpDefn[1]__2____28________ = 1 ->
   Assign['RVY']____ <= 0
 MinInGrpDefn[4]__2____26____4___: MinInGrpDefn[4]__2____26________ = 1 ->
   Assign['RCR']____ <= 3
 MinInGrpDefn[3]__2____26____4___: MinInGrpDefn[3]__2____26________ = 1 ->
   Assign['RCR']____ <= 2
 MinInGrpDefn[5]__2____26____4___: MinInGrpDefn[5]__2____26________ = 1 ->
   Assign['RCR']____ <= 4
 MinInGrpDefn[2]__2____26____4___: MinInGrpDefn[2]__2____26________ = 1 ->
   Assign['RCR']____ <= 1
 MinInGrpDefn[1]__2____26____4___: MinInGrpDefn[1]__2____26________ = 1 ->
   Assign['RCR']____ <= 0
 MinInGrpDefn[4]__2____24____4___: MinInGrpDefn[4]__2____24________ = 1 ->
   Assign['SCR']____ <= 3
 MinInGrpDefn[3]__2____24____4___: MinInGrpDefn[3]__2____24________ = 1 ->
   Assign['SCR']____ <= 2
 MinInGrpDefn[5]__2____24____4___: MinInGrpDefn[5]__2____24________ = 1 ->
   Assign['SCR']____ <= 4
 MinInGrpDefn[2]__2____24____4___: MinInGrpDefn[2]__2____24________ = 1 ->
   Assign['SCR']____ <= 1
 MinInGrpDefn[1]__2____24____4___: MinInGrpDefn[1]__2____24________ = 1 ->
   Assign['SCR']____ <= 0
 MinInGrpDefn[4]__2____22____4___: MinInGrpDefn[4]__2____22________ = 1 ->
   Assign['MBT']____ <= 3
 MinInGrpDefn[3]__2____22____4___: MinInGrpDefn[3]__2____22________ = 1 ->
   Assign['MBT']____ <= 2
 MinInGrpDefn[5]__2____22____4___: MinInGrpDefn[5]__2____22________ = 1 ->
   Assign['MBT']____ <= 4
 MinInGrpDefn[2]__2____22____4___: MinInGrpDefn[2]__2____22________ = 1 ->
   Assign['MBT']____ <= 1
 MinInGrpDefn[1]__2____22____4___: MinInGrpDefn[1]__2____22________ = 1 ->
   Assign['MBT']____ <= 0
 MinInGrpDefn[4]__2____20____4___: MinInGrpDefn[4]__2____20________ = 1 ->
   Assign['JWR']____ <= 3
 MinInGrpDefn[3]__2____20____4___: MinInGrpDefn[3]__2____20________ = 1 ->
   Assign['JWR']____ <= 2
 MinInGrpDefn[5]__2____20____4___: MinInGrpDefn[5]__2____20________ = 1 ->
   Assign['JWR']____ <= 4
 MinInGrpDefn[2]__2____20____4___: MinInGrpDefn[2]__2____20________ = 1 ->
   Assign['JWR']____ <= 1
 MinInGrpDefn[1]__2____20____4___: MinInGrpDefn[1]__2____20________ = 1 ->
   Assign['JWR']____ <= 0
 MinInGrpDefn[4]__2____18____4___: MinInGrpDefn[4]__2____18________ = 1 ->
   Assign['CWU']____ <= 3
 MinInGrpDefn[3]__2____18____4___: MinInGrpDefn[3]__2____18________ = 1 ->
   Assign['CWU']____ <= 2
 MinInGrpDefn[5]__2____18____4___: MinInGrpDefn[5]__2____18________ = 1 ->
   Assign['CWU']____ <= 4
 MinInGrpDefn[2]__2____18____4___: MinInGrpDefn[2]__2____18________ = 1 ->
   Assign['CWU']____ <= 1
 MinInGrpDefn[1]__2____18____4___: MinInGrpDefn[1]__2____18________ = 1 ->
   Assign['CWU']____ <= 0
 MinInGrpDefn[4]__2____16____4___: MinInGrpDefn[4]__2____16________ = 1 ->
   Assign['ACI']____ <= 3
 MinInGrpDefn[3]__2____16____4___: MinInGrpDefn[3]__2____16________ = 1 ->
   Assign['ACI']____ <= 2
 MinInGrpDefn[5]__2____16____4___: MinInGrpDefn[5]__2____16________ = 1 ->
   Assign['ACI']____ <= 4
 MinInGrpDefn[2]__2____16____4___: MinInGrpDefn[2]__2____16________ = 1 ->
   Assign['ACI']____ <= 1
 MinInGrpDefn[1]__2____16____4___: MinInGrpDefn[1]__2____16________ = 1 ->
   Assign['ACI']____ <= 0
 MinInGrpDefn[4]__2____14____4___: MinInGrpDefn[4]__2____14________ = 1 ->
   Assign['BLE']____ <= 3
 MinInGrpDefn[3]__2____14____4___: MinInGrpDefn[3]__2____14________ = 1 ->
   Assign['BLE']____ <= 2
 MinInGrpDefn[5]__2____14____4___: MinInGrpDefn[5]__2____14________ = 1 ->
   Assign['BLE']____ <= 4
 MinInGrpDefn[2]__2____14____4___: MinInGrpDefn[2]__2____14________ = 1 ->
   Assign['BLE']____ <= 1
 MinInGrpDefn[1]__2____14____4___: MinInGrpDefn[1]__2____14________ = 1 ->
   Assign['BLE']____ <= 0
 MinInGrpDefn[4]__2____12____4___: MinInGrpDefn[4]__2____12________ = 1 ->
   Assign['SDN']____ <= 3
 MinInGrpDefn[3]__2____12____4___: MinInGrpDefn[3]__2____12________ = 1 ->
   Assign['SDN']____ <= 2
 MinInGrpDefn[5]__2____12____4___: MinInGrpDefn[5]__2____12________ = 1 ->
   Assign['SDN']____ <= 4
 MinInGrpDefn[2]__2____12____4___: MinInGrpDefn[2]__2____12________ = 1 ->
   Assign['SDN']____ <= 1
 MinInGrpDefn[1]__2____12____4___: MinInGrpDefn[1]__2____12________ = 1 ->
   Assign['SDN']____ <= 0
 MinInGrpDefn[4]__2____10____4___: MinInGrpDefn[4]__2____10________ = 1 ->
   Assign['JLN']____ <= 3
 MinInGrpDefn[3]__2____10____4___: MinInGrpDefn[3]__2____10________ = 1 ->
   Assign['JLN']____ <= 2
 MinInGrpDefn[5]__2____10____4___: MinInGrpDefn[5]__2____10________ = 1 ->
   Assign['JLN']____ <= 4
 MinInGrpDefn[2]__2____10____4___: MinInGrpDefn[2]__2____10________ = 1 ->
   Assign['JLN']____ <= 1
 MinInGrpDefn[1]__2____10____4___: MinInGrpDefn[1]__2____10________ = 1 ->
   Assign['JLN']____ <= 0
 MinInGrpDefn[4]__2____6____4___: MinInGrpDefn[4]__2____6________ = 1 ->
   Assign['JAE']____ <= 3
 MinInGrpDefn[3]__2____6____4___: MinInGrpDefn[3]__2____6________ = 1 ->
   Assign['JAE']____ <= 2
 MinInGrpDefn[5]__2____6____4___: MinInGrpDefn[5]__2____6________ = 1 ->
   Assign['JAE']____ <= 4
 MinInGrpDefn[2]__2____6____4___: MinInGrpDefn[2]__2____6________ = 1 ->
   Assign['JAE']____ <= 1
 MinInGrpDefn[1]__2____6____4___: MinInGrpDefn[1]__2____6________ = 1 ->
   Assign['JAE']____ <= 0
 MinInGrpDefn[4]__2____2____4___: MinInGrpDefn[4]__2____2________ = 1 ->
   Assign['SML']____ <= 3
 MinInGrpDefn[3]__2____2____4___: MinInGrpDefn[3]__2____2________ = 1 ->
   Assign['SML']____ <= 2
 MinInGrpDefn[5]__2____2____4___: MinInGrpDefn[5]__2____2________ = 1 ->
   Assign['SML']____ <= 4
 MinInGrpDefn[2]__2____2____4___: MinInGrpDefn[2]__2____2________ = 1 ->
   Assign['SML']____ <= 1
 MinInGrpDefn[1]__2____2____4___: MinInGrpDefn[1]__2____2________ = 1 ->
   Assign['SML']____ <= 0
 MinInGrpDefn[4]__2____8____4___: MinInGrpDefn[4]__2____8________ = 1 ->
   Assign['FTR']____ <= 3
 MinInGrpDefn[3]__2____8____4___: MinInGrpDefn[3]__2____8________ = 1 ->
   Assign['FTR']____ <= 2
 MinInGrpDefn[5]__2____8____4___: MinInGrpDefn[5]__2____8________ = 1 ->
   Assign['FTR']____ <= 4
 MinInGrpDefn[2]__2____8____4___: MinInGrpDefn[2]__2____8________ = 1 ->
   Assign['FTR']____ <= 1
 MinInGrpDefn[1]__2____8____4___: MinInGrpDefn[1]__2____8________ = 1 ->
   Assign['FTR']____ <= 0
 MinInGrpDefn[4]__2____50____5___: MinInGrpDefn[4]__2____50____2_____ = 1 ->
   Assign['JEN']____ >= 5
 MinInGrpDefn[3]__2____50____5___: MinInGrpDefn[3]__2____50____2_____ = 1 ->
   Assign['JEN']____ >= 4
 MinInGrpDefn[5]__2____50____5___: MinInGrpDefn[5]__2____50____2_____ = 1 ->
   Assign['JEN']____ >= 6
 MinInGrpDefn[2]__2____50____5___: MinInGrpDefn[2]__2____50____2_____ = 1 ->
   Assign['JEN']____ >= 3
 MinInGrpDefn[1]__2____50____5___: MinInGrpDefn[1]__2____50____2_____ = 1 ->
   Assign['JEN']____ >= 2
 MinInGrpDefn[4]__2____48____5___: MinInGrpDefn[4]__2____48____2_____ = 1 ->
   Assign['PSN']____ >= 5
 MinInGrpDefn[3]__2____48____5___: MinInGrpDefn[3]__2____48____2_____ = 1 ->
   Assign['PSN']____ >= 4
 MinInGrpDefn[5]__2____48____5___: MinInGrpDefn[5]__2____48____2_____ = 1 ->
   Assign['PSN']____ >= 6
 MinInGrpDefn[2]__2____48____5___: MinInGrpDefn[2]__2____48____2_____ = 1 ->
   Assign['PSN']____ >= 3
 MinInGrpDefn[1]__2____48____5___: MinInGrpDefn[1]__2____48____2_____ = 1 ->
   Assign['PSN']____ >= 2
 MinInGrpDefn[4]__2____46____5___: MinInGrpDefn[4]__2____46____2_____ = 1 ->
   Assign['PSO']____ >= 5
 MinInGrpDefn[3]__2____46____5___: MinInGrpDefn[3]__2____46____2_____ = 1 ->
   Assign['PSO']____ >= 4
 MinInGrpDefn[5]__2____46____5___: MinInGrpDefn[5]__2____46____2_____ = 1 ->
   Assign['PSO']____ >= 6
 MinInGrpDefn[2]__2____46____5___: MinInGrpDefn[2]__2____46____2_____ = 1 ->
   Assign['PSO']____ >= 3
 MinInGrpDefn[1]__2____46____5___: MinInGrpDefn[1]__2____46____2_____ = 1 ->
   Assign['PSO']____ >= 2
 MinInGrpDefn[4]__2____44____5___: MinInGrpDefn[4]__2____44____2_____ = 1 ->
   Assign['JMA']____ >= 5
 MinInGrpDefn[3]__2____44____5___: MinInGrpDefn[3]__2____44____2_____ = 1 ->
   Assign['JMA']____ >= 4
 MinInGrpDefn[5]__2____44____5___: MinInGrpDefn[5]__2____44____2_____ = 1 ->
   Assign['JMA']____ >= 6
 MinInGrpDefn[2]__2____44____5___: MinInGrpDefn[2]__2____44____2_____ = 1 ->
   Assign['JMA']____ >= 3
 MinInGrpDefn[1]__2____44____5___: MinInGrpDefn[1]__2____44____2_____ = 1 ->
   Assign['JMA']____ >= 2
 MinInGrpDefn[4]__2____42____5___: MinInGrpDefn[4]__2____42____2_____ = 1 ->
   Assign['PSA']____ >= 5
 MinInGrpDefn[3]__2____42____5___: MinInGrpDefn[3]__2____42____2_____ = 1 ->
   Assign['PSA']____ >= 4
 MinInGrpDefn[5]__2____42____5___: MinInGrpDefn[5]__2____42____2_____ = 1 ->
   Assign['PSA']____ >= 6
 MinInGrpDefn[2]__2____42____5___: MinInGrpDefn[2]__2____42____2_____ = 1 ->
   Assign['PSA']____ >= 3
 MinInGrpDefn[1]__2____42____5___: MinInGrpDefn[1]__2____42____2_____ = 1 ->
   Assign['PSA']____ >= 2
 MinInGrpDefn[4]__2____40____5___: MinInGrpDefn[4]__2____40____2_____ = 1 ->
   Assign['RGE']____ >= 5
 MinInGrpDefn[3]__2____40____5___: MinInGrpDefn[3]__2____40____2_____ = 1 ->
   Assign['RGE']____ >= 4
 MinInGrpDefn[5]__2____40____5___: MinInGrpDefn[5]__2____40____2_____ = 1 ->
   Assign['RGE']____ >= 6
 MinInGrpDefn[2]__2____40____5___: MinInGrpDefn[2]__2____40____2_____ = 1 ->
   Assign['RGE']____ >= 3
 MinInGrpDefn[1]__2____40____5___: MinInGrpDefn[1]__2____40____2_____ = 1 ->
   Assign['RGE']____ >= 2
 MinInGrpDefn[4]__2____38____5___: MinInGrpDefn[4]__2____38____2_____ = 1 ->
   Assign['FEE']____ >= 5
 MinInGrpDefn[3]__2____38____5___: MinInGrpDefn[3]__2____38____2_____ = 1 ->
   Assign['FEE']____ >= 4
 MinInGrpDefn[5]__2____38____5___: MinInGrpDefn[5]__2____38____2_____ = 1 ->
   Assign['FEE']____ >= 6
 MinInGrpDefn[2]__2____38____5___: MinInGrpDefn[2]__2____38____2_____ = 1 ->
   Assign['FEE']____ >= 3
 MinInGrpDefn[1]__2____38____5___: MinInGrpDefn[1]__2____38____2_____ = 1 ->
   Assign['FEE']____ >= 2
 MinInGrpDefn[4]__2____36____5___: MinInGrpDefn[4]__2____36____2_____ = 1 ->
   Assign['RHA']____ >= 5
 MinInGrpDefn[3]__2____36____5___: MinInGrpDefn[3]__2____36____2_____ = 1 ->
   Assign['RHA']____ >= 4
 MinInGrpDefn[5]__2____36____5___: MinInGrpDefn[5]__2____36____2_____ = 1 ->
   Assign['RHA']____ >= 6
 MinInGrpDefn[2]__2____36____5___: MinInGrpDefn[2]__2____36____2_____ = 1 ->
   Assign['RHA']____ >= 3
 MinInGrpDefn[1]__2____36____5___: MinInGrpDefn[1]__2____36____2_____ = 1 ->
   Assign['RHA']____ >= 2
 MinInGrpDefn[4]__2____34____5___: MinInGrpDefn[4]__2____34____2_____ = 1 ->
   Assign['VMA']____ >= 5
 MinInGrpDefn[3]__2____34____5___: MinInGrpDefn[3]__2____34____2_____ = 1 ->
   Assign['VMA']____ >= 4
 MinInGrpDefn[5]__2____34____5___: MinInGrpDefn[5]__2____34____2_____ = 1 ->
   Assign['VMA']____ >= 6
 MinInGrpDefn[2]__2____34____5___: MinInGrpDefn[2]__2____34____2_____ = 1 ->
   Assign['VMA']____ >= 3
 MinInGrpDefn[1]__2____34____5___: MinInGrpDefn[1]__2____34____2_____ = 1 ->
   Assign['VMA']____ >= 2
 MinInGrpDefn[4]__2____32____5___: MinInGrpDefn[4]__2____32____2_____ = 1 ->
   Assign['PST']____ >= 5
 MinInGrpDefn[3]__2____32____5___: MinInGrpDefn[3]__2____32____2_____ = 1 ->
   Assign['PST']____ >= 4
 MinInGrpDefn[5]__2____32____5___: MinInGrpDefn[5]__2____32____2_____ = 1 ->
   Assign['PST']____ >= 6
 MinInGrpDefn[2]__2____32____5___: MinInGrpDefn[2]__2____32____2_____ = 1 ->
   Assign['PST']____ >= 3
 MinInGrpDefn[1]__2____32____5___: MinInGrpDefn[1]__2____32____2_____ = 1 ->
   Assign['PST']____ >= 2
 MinInGrpDefn[4]__2____28____5___: MinInGrpDefn[4]__2____28____2_____ = 1 ->
   Assign['RVY']____ >= 5
 MinInGrpDefn[3]__2____28____5___: MinInGrpDefn[3]__2____28____2_____ = 1 ->
   Assign['RVY']____ >= 4
 MinInGrpDefn[5]__2____28____5___: MinInGrpDefn[5]__2____28____2_____ = 1 ->
   Assign['RVY']____ >= 6
 MinInGrpDefn[2]__2____28____5___: MinInGrpDefn[2]__2____28____2_____ = 1 ->
   Assign['RVY']____ >= 3
 MinInGrpDefn[1]__2____28____5___: MinInGrpDefn[1]__2____28____2_____ = 1 ->
   Assign['RVY']____ >= 2
 MinInGrpDefn[4]__2____26____5___: MinInGrpDefn[4]__2____26____2_____ = 1 ->
   Assign['RCR']____ >= 5
 MinInGrpDefn[3]__2____26____5___: MinInGrpDefn[3]__2____26____2_____ = 1 ->
   Assign['RCR']____ >= 4
 MinInGrpDefn[5]__2____26____5___: MinInGrpDefn[5]__2____26____2_____ = 1 ->
   Assign['RCR']____ >= 6
 MinInGrpDefn[2]__2____26____5___: MinInGrpDefn[2]__2____26____2_____ = 1 ->
   Assign['RCR']____ >= 3
 MinInGrpDefn[1]__2____26____5___: MinInGrpDefn[1]__2____26____2_____ = 1 ->
   Assign['RCR']____ >= 2
 MinInGrpDefn[4]__2____24____5___: MinInGrpDefn[4]__2____24____2_____ = 1 ->
   Assign['SCR']____ >= 5
 MinInGrpDefn[3]__2____24____5___: MinInGrpDefn[3]__2____24____2_____ = 1 ->
   Assign['SCR']____ >= 4
 MinInGrpDefn[5]__2____24____5___: MinInGrpDefn[5]__2____24____2_____ = 1 ->
   Assign['SCR']____ >= 6
 MinInGrpDefn[2]__2____24____5___: MinInGrpDefn[2]__2____24____2_____ = 1 ->
   Assign['SCR']____ >= 3
 MinInGrpDefn[1]__2____24____5___: MinInGrpDefn[1]__2____24____2_____ = 1 ->
   Assign['SCR']____ >= 2
 MinInGrpDefn[4]__2____22____5___: MinInGrpDefn[4]__2____22____2_____ = 1 ->
   Assign['MBT']____ >= 5
 MinInGrpDefn[3]__2____22____5___: MinInGrpDefn[3]__2____22____2_____ = 1 ->
   Assign['MBT']____ >= 4
 MinInGrpDefn[5]__2____22____5___: MinInGrpDefn[5]__2____22____2_____ = 1 ->
   Assign['MBT']____ >= 6
 MinInGrpDefn[2]__2____22____5___: MinInGrpDefn[2]__2____22____2_____ = 1 ->
   Assign['MBT']____ >= 3
 MinInGrpDefn[1]__2____22____5___: MinInGrpDefn[1]__2____22____2_____ = 1 ->
   Assign['MBT']____ >= 2
 MinInGrpDefn[4]__2____20____5___: MinInGrpDefn[4]__2____20____2_____ = 1 ->
   Assign['JWR']____ >= 5
 MinInGrpDefn[3]__2____20____5___: MinInGrpDefn[3]__2____20____2_____ = 1 ->
   Assign['JWR']____ >= 4
 MinInGrpDefn[5]__2____20____5___: MinInGrpDefn[5]__2____20____2_____ = 1 ->
   Assign['JWR']____ >= 6
 MinInGrpDefn[2]__2____20____5___: MinInGrpDefn[2]__2____20____2_____ = 1 ->
   Assign['JWR']____ >= 3
 MinInGrpDefn[1]__2____20____5___: MinInGrpDefn[1]__2____20____2_____ = 1 ->
   Assign['JWR']____ >= 2
 MinInGrpDefn[4]__2____18____5___: MinInGrpDefn[4]__2____18____2_____ = 1 ->
   Assign['CWU']____ >= 5
 MinInGrpDefn[3]__2____18____5___: MinInGrpDefn[3]__2____18____2_____ = 1 ->
   Assign['CWU']____ >= 4
 MinInGrpDefn[5]__2____18____5___: MinInGrpDefn[5]__2____18____2_____ = 1 ->
   Assign['CWU']____ >= 6
 MinInGrpDefn[2]__2____18____5___: MinInGrpDefn[2]__2____18____2_____ = 1 ->
   Assign['CWU']____ >= 3
 MinInGrpDefn[1]__2____18____5___: MinInGrpDefn[1]__2____18____2_____ = 1 ->
   Assign['CWU']____ >= 2
 MinInGrpDefn[4]__2____16____5___: MinInGrpDefn[4]__2____16____2_____ = 1 ->
   Assign['ACI']____ >= 5
 MinInGrpDefn[3]__2____16____5___: MinInGrpDefn[3]__2____16____2_____ = 1 ->
   Assign['ACI']____ >= 4
 MinInGrpDefn[5]__2____16____5___: MinInGrpDefn[5]__2____16____2_____ = 1 ->
   Assign['ACI']____ >= 6
 MinInGrpDefn[2]__2____16____5___: MinInGrpDefn[2]__2____16____2_____ = 1 ->
   Assign['ACI']____ >= 3
 MinInGrpDefn[1]__2____16____5___: MinInGrpDefn[1]__2____16____2_____ = 1 ->
   Assign['ACI']____ >= 2
 MinInGrpDefn[4]__2____14____5___: MinInGrpDefn[4]__2____14____2_____ = 1 ->
   Assign['BLE']____ >= 5
 MinInGrpDefn[3]__2____14____5___: MinInGrpDefn[3]__2____14____2_____ = 1 ->
   Assign['BLE']____ >= 4
 MinInGrpDefn[5]__2____14____5___: MinInGrpDefn[5]__2____14____2_____ = 1 ->
   Assign['BLE']____ >= 6
 MinInGrpDefn[2]__2____14____5___: MinInGrpDefn[2]__2____14____2_____ = 1 ->
   Assign['BLE']____ >= 3
 MinInGrpDefn[1]__2____14____5___: MinInGrpDefn[1]__2____14____2_____ = 1 ->
   Assign['BLE']____ >= 2
 MinInGrpDefn[4]__2____12____5___: MinInGrpDefn[4]__2____12____2_____ = 1 ->
   Assign['SDN']____ >= 5
 MinInGrpDefn[3]__2____12____5___: MinInGrpDefn[3]__2____12____2_____ = 1 ->
   Assign['SDN']____ >= 4
 MinInGrpDefn[5]__2____12____5___: MinInGrpDefn[5]__2____12____2_____ = 1 ->
   Assign['SDN']____ >= 6
 MinInGrpDefn[2]__2____12____5___: MinInGrpDefn[2]__2____12____2_____ = 1 ->
   Assign['SDN']____ >= 3
 MinInGrpDefn[1]__2____12____5___: MinInGrpDefn[1]__2____12____2_____ = 1 ->
   Assign['SDN']____ >= 2
 MinInGrpDefn[4]__2____10____5___: MinInGrpDefn[4]__2____10____2_____ = 1 ->
   Assign['JLN']____ >= 5
 MinInGrpDefn[3]__2____10____5___: MinInGrpDefn[3]__2____10____2_____ = 1 ->
   Assign['JLN']____ >= 4
 MinInGrpDefn[5]__2____10____5___: MinInGrpDefn[5]__2____10____2_____ = 1 ->
   Assign['JLN']____ >= 6
 MinInGrpDefn[2]__2____10____5___: MinInGrpDefn[2]__2____10____2_____ = 1 ->
   Assign['JLN']____ >= 3
 MinInGrpDefn[1]__2____10____5___: MinInGrpDefn[1]__2____10____2_____ = 1 ->
   Assign['JLN']____ >= 2
 MinInGrpDefn[4]__2____6____5___: MinInGrpDefn[4]__2____6____2_____ = 1 ->
   Assign['JAE']____ >= 5
 MinInGrpDefn[3]__2____6____5___: MinInGrpDefn[3]__2____6____2_____ = 1 ->
   Assign['JAE']____ >= 4
 MinInGrpDefn[5]__2____6____5___: MinInGrpDefn[5]__2____6____2_____ = 1 ->
   Assign['JAE']____ >= 6
 MinInGrpDefn[2]__2____6____5___: MinInGrpDefn[2]__2____6____2_____ = 1 ->
   Assign['JAE']____ >= 3
 MinInGrpDefn[1]__2____6____5___: MinInGrpDefn[1]__2____6____2_____ = 1 ->
   Assign['JAE']____ >= 2
 MinInGrpDefn[4]__2____2____5___: MinInGrpDefn[4]__2____2____2_____ = 1 ->
   Assign['SML']____ >= 5
 MinInGrpDefn[3]__2____2____5___: MinInGrpDefn[3]__2____2____2_____ = 1 ->
   Assign['SML']____ >= 4
 MinInGrpDefn[5]__2____2____5___: MinInGrpDefn[5]__2____2____2_____ = 1 ->
   Assign['SML']____ >= 6
 MinInGrpDefn[2]__2____2____5___: MinInGrpDefn[2]__2____2____2_____ = 1 ->
   Assign['SML']____ >= 3
 MinInGrpDefn[1]__2____2____5___: MinInGrpDefn[1]__2____2____2_____ = 1 ->
   Assign['SML']____ >= 2
 MinInGrpDefn[4]__2____8____5___: MinInGrpDefn[4]__2____8____2_____ = 1 ->
   Assign['FTR']____ >= 5
 MinInGrpDefn[3]__2____8____5___: MinInGrpDefn[3]__2____8____2_____ = 1 ->
   Assign['FTR']____ >= 4
 MinInGrpDefn[5]__2____8____5___: MinInGrpDefn[5]__2____8____2_____ = 1 ->
   Assign['FTR']____ >= 6
 MinInGrpDefn[2]__2____8____5___: MinInGrpDefn[2]__2____8____2_____ = 1 ->
   Assign['FTR']____ >= 3
 MinInGrpDefn[1]__2____8____5___: MinInGrpDefn[1]__2____8____2_____ = 1 ->
   Assign['FTR']____ >= 2
Bounds
 1 <= Assign['SML']____ <= 5
 1 <= Assign['BCG']____ <= 5
 1 <= Assign['JAE']____ <= 5
 1 <= Assign['FTR']____ <= 5
 1 <= Assign['JLN']____ <= 5
 1 <= Assign['SDN']____ <= 5
 1 <= Assign['BLE']____ <= 5
 1 <= Assign['ACI']____ <= 5
 1 <= Assign['CWU']____ <= 5
 1 <= Assign['JWR']____ <= 5
 1 <= Assign['MBT']____ <= 5
 1 <= Assign['SCR']____ <= 5
 1 <= Assign['RCR']____ <= 5
 1 <= Assign['RVY']____ <= 5
 1 <= Assign['KWN']____ <= 5
 1 <= Assign['PST']____ <= 5
 1 <= Assign['VMA']____ <= 5
 1 <= Assign['RHA']____ <= 5
 1 <= Assign['FEE']____ <= 5
 1 <= Assign['RGE']____ <= 5
 1 <= Assign['PSA']____ <= 5
 1 <= Assign['JMA']____ <= 5
 1 <= Assign['PSO']____ <= 5
 1 <= Assign['PSN']____ <= 5
 1 <= Assign['JEN']____ <= 5
 1 <= Assign['RSN']____ <= 5
 -infinity <= MinInGrp____ <= 5
 MaxInGrp____ >= 6
 -infinity <= MinType['dept','NNE']____ <= 4
 -infinity <= MinType['dept','NNW']____ <= 0
 -infinity <= MinType['dept','WSW']____ <= 0
 -infinity <= MinType['dept','SSE']____ <= 0
 -infinity <= MinType['loc','Springfield']____ <= 1
 -infinity <= MinType['loc','Urbana']____ <= 0
 -infinity <= MinType['loc','Peoria']____ <= 2
 -infinity <= MinType['loc','Carbondale']____ <= 0
 -infinity <= MinType['loc','Joliet']____ <= 0
 -infinity <= MinType['rate','A']____ <= 4
 -infinity <= MinType['rate','B']____ <= 0
 -infinity <= MinType['title','Assistant']____ <= 2
 -infinity <= MinType['title','Adjunct']____ <= 2
 -infinity <= MinType['title','Deputy']____ <= 0
 MaxType['dept','NNE']____ >= 4
 MaxType['dept','NNW']____ >= 1
 MaxType['dept','WSW']____ >= 1
 MaxType['dept','SSE']____ >= 1
 MaxType['loc','Springfield']____ >= 2
 MaxType['loc','Urbana']____ >= 1
 MaxType['loc','Peoria']____ >= 3
 MaxType['loc','Carbondale']____ >= 1
 MaxType['loc','Joliet']____ >= 1
 MaxType['rate','A']____ >= 5
 MaxType['rate','B']____ >= 1
 MaxType['title','Assistant']____ >= 3
 MaxType['title','Adjunct']____ >= 3
 MaxType['title','Deputy']____ >= 1
 MinInGrpDefn[1]______ <= 26
 MinInGrpDefn[2]______ <= 26
 MinInGrpDefn[3]______ <= 26
 MinInGrpDefn[4]______ <= 26
 MinInGrpDefn[5]______ <= 26
 MinTypeDefn[1,'dept','NNE']______ <= 20
 MinTypeDefn[1,'dept','NNW']______ <= 3
 MinTypeDefn[1,'dept','WSW']______ <= 1
 MinTypeDefn[1,'dept','WSW']__3_____ = 1
 MinTypeDefn[1,'dept','WSW']__4_____ = 0
 MinTypeDefn[1,'dept','WSW']__5_____ = 0
 MinTypeDefn[1,'dept','SSE']______ <= 2
 MinTypeDefn[1,'loc','Springfield']______ <= 8
 MinTypeDefn[1,'loc','Urbana']______ <= 1
 MinTypeDefn[1,'loc','Urbana']__5_____ = 0
 MinTypeDefn[1,'loc','Peoria']______ <= 14
 MinTypeDefn[1,'loc','Carbondale']______ <= 2
 MinTypeDefn[1,'loc','Joliet']______ <= 1
 MinTypeDefn[1,'loc','Joliet']__5_____ = 0
 MinTypeDefn[1,'rate','A']______ <= 24
 MinTypeDefn[1,'rate','B']______ <= 2
 MinTypeDefn[1,'title','Assistant']______ <= 12
 MinTypeDefn[1,'title','Adjunct']______ <= 11
 MinTypeDefn[1,'title','Deputy']______ <= 3
 MinTypeDefn[2,'dept','NNE']______ <= 20
 MinTypeDefn[2,'dept','NNW']______ <= 3
 MinTypeDefn[2,'dept','WSW']______ <= 1
 MinTypeDefn[2,'dept','WSW']__5_____ = 0
 MinTypeDefn[2,'dept','SSE']______ <= 2
 MinTypeDefn[2,'loc','Springfield']______ <= 8
 MinTypeDefn[2,'loc','Urbana']______ <= 1
 MinTypeDefn[2,'loc','Urbana']__5_____ = 0
 MinTypeDefn[2,'loc','Peoria']______ <= 14
 MinTypeDefn[2,'loc','Carbondale']______ <= 2
 MinTypeDefn[2,'loc','Joliet']______ <= 1
 MinTypeDefn[2,'loc','Joliet']__5_____ = 0
 MinTypeDefn[2,'rate','A']______ <= 24
 MinTypeDefn[2,'rate','B']______ <= 2
 MinTypeDefn[2,'title','Assistant']______ <= 12
 MinTypeDefn[2,'title','Adjunct']______ <= 11
 MinTypeDefn[2,'title','Deputy']______ <= 3
 MinTypeDefn[3,'dept','NNE']______ <= 20
 MinTypeDefn[3,'dept','NNW']______ <= 3
 MinTypeDefn[3,'dept','WSW']______ <= 1
 MinTypeDefn[3,'dept','WSW']__5_____ = 0
 MinTypeDefn[3,'dept','SSE']______ <= 2
 MinTypeDefn[3,'loc','Springfield']______ <= 8
 MinTypeDefn[3,'loc','Urbana']______ <= 1
 MinTypeDefn[3,'loc','Urbana']__5_____ = 0
 MinTypeDefn[3,'loc','Peoria']______ <= 14
 MinTypeDefn[3,'loc','Carbondale']______ <= 2
 MinTypeDefn[3,'loc','Joliet']______ <= 1
 MinTypeDefn[3,'loc','Joliet']__5_____ = 0
 MinTypeDefn[3,'rate','A']______ <= 24
 MinTypeDefn[3,'rate','B']______ <= 2
 MinTypeDefn[3,'title','Assistant']______ <= 12
 MinTypeDefn[3,'title','Adjunct']______ <= 11
 MinTypeDefn[3,'title','Deputy']______ <= 3
 MinTypeDefn[4,'dept','NNE']______ <= 20
 MinTypeDefn[4,'dept','NNW']______ <= 3
 MinTypeDefn[4,'dept','WSW']______ <= 1
 MinTypeDefn[4,'dept','WSW']__5_____ = 0
 MinTypeDefn[4,'dept','SSE']______ <= 2
 MinTypeDefn[4,'loc','Springfield']______ <= 8
 MinTypeDefn[4,'loc','Urbana']______ <= 1
 MinTypeDefn[4,'loc','Urbana']__5_____ = 0
 MinTypeDefn[4,'loc','Peoria']______ <= 14
 MinTypeDefn[4,'loc','Carbondale']______ <= 2
 MinTypeDefn[4,'loc','Joliet']______ <= 1
 MinTypeDefn[4,'loc','Joliet']__5_____ = 0
 MinTypeDefn[4,'rate','A']______ <= 24
 MinTypeDefn[4,'rate','B']______ <= 2
 MinTypeDefn[4,'title','Assistant']______ <= 12
 MinTypeDefn[4,'title','Adjunct']______ <= 11
 MinTypeDefn[4,'title','Deputy']______ <= 3
 MinTypeDefn[5,'dept','NNE']______ <= 20
 MinTypeDefn[5,'dept','NNW']______ <= 3
 MinTypeDefn[5,'dept','WSW']______ <= 1
 MinTypeDefn[5,'dept','WSW']__5_____ = 0
 MinTypeDefn[5,'dept','SSE']______ <= 2
 MinTypeDefn[5,'loc','Springfield']______ <= 8
 MinTypeDefn[5,'loc','Urbana']______ <= 1
 MinTypeDefn[5,'loc','Urbana']__5_____ = 0
 MinTypeDefn[5,'loc','Peoria']______ <= 14
 MinTypeDefn[5,'loc','Carbondale']______ <= 2
 MinTypeDefn[5,'loc','Joliet']______ <= 1
 MinTypeDefn[5,'loc','Joliet']__5_____ = 0
 MinTypeDefn[5,'rate','A']______ <= 24
 MinTypeDefn[5,'rate','B']______ <= 2
 MinTypeDefn[5,'title','Assistant']______ <= 12
 MinTypeDefn[5,'title','Adjunct']______ <= 11
 MinTypeDefn[5,'title','Deputy']______ <= 3
 MinInGrpDefn[1]__2________ <= 1
 MinInGrpDefn[1]__2____5_____ <= 1
 MinInGrpDefn[1]__2____7_____ <= 1
 MinInGrpDefn[1]__2____9_____ <= 1
 MinInGrpDefn[1]__2____11_____ <= 1
 MinInGrpDefn[1]__2____13_____ <= 1
 MinInGrpDefn[1]__2____15_____ <= 1
 MinInGrpDefn[1]__2____17_____ <= 1
 MinInGrpDefn[1]__2____19_____ <= 1
 MinInGrpDefn[1]__2____21_____ <= 1
 MinInGrpDefn[1]__2____23_____ <= 1
 MinInGrpDefn[1]__2____25_____ <= 1
 MinInGrpDefn[1]__2____27_____ <= 1
 MinInGrpDefn[1]__2____31_____ <= 1
 MinInGrpDefn[1]__2____33_____ <= 1
 MinInGrpDefn[1]__2____35_____ <= 1
 MinInGrpDefn[1]__2____37_____ <= 1
 MinInGrpDefn[1]__2____39_____ <= 1
 MinInGrpDefn[1]__2____41_____ <= 1
 MinInGrpDefn[1]__2____43_____ <= 1
 MinInGrpDefn[1]__2____45_____ <= 1
 MinInGrpDefn[1]__2____47_____ <= 1
 MinInGrpDefn[1]__2____49_____ <= 1
 MinInGrpDefn[2]__2________ <= 1
 MinInGrpDefn[2]__2____5_____ <= 1
 MinInGrpDefn[2]__2____7_____ <= 1
 MinInGrpDefn[2]__2____9_____ <= 1
 MinInGrpDefn[2]__2____11_____ <= 1
 MinInGrpDefn[2]__2____13_____ <= 1
 MinInGrpDefn[2]__2____15_____ <= 1
 MinInGrpDefn[2]__2____17_____ <= 1
 MinInGrpDefn[2]__2____19_____ <= 1
 MinInGrpDefn[2]__2____21_____ <= 1
 MinInGrpDefn[2]__2____23_____ <= 1
 MinInGrpDefn[2]__2____25_____ <= 1
 MinInGrpDefn[2]__2____27_____ <= 1
 MinInGrpDefn[2]__2____31_____ <= 1
 MinInGrpDefn[2]__2____33_____ <= 1
 MinInGrpDefn[2]__2____35_____ <= 1
 MinInGrpDefn[2]__2____37_____ <= 1
 MinInGrpDefn[2]__2____39_____ <= 1
 MinInGrpDefn[2]__2____41_____ <= 1
 MinInGrpDefn[2]__2____43_____ <= 1
 MinInGrpDefn[2]__2____45_____ <= 1
 MinInGrpDefn[2]__2____47_____ <= 1
 MinInGrpDefn[2]__2____49_____ <= 1
 MinInGrpDefn[3]__2________ <= 1
 MinInGrpDefn[3]__2____5_____ <= 1
 MinInGrpDefn[3]__2____7_____ <= 1
 MinInGrpDefn[3]__2____9_____ <= 1
 MinInGrpDefn[3]__2____11_____ <= 1
 MinInGrpDefn[3]__2____13_____ <= 1
 MinInGrpDefn[3]__2____15_____ <= 1
 MinInGrpDefn[3]__2____17_____ <= 1
 MinInGrpDefn[3]__2____19_____ <= 1
 MinInGrpDefn[3]__2____21_____ <= 1
 MinInGrpDefn[3]__2____23_____ <= 1
 MinInGrpDefn[3]__2____25_____ <= 1
 MinInGrpDefn[3]__2____27_____ <= 1
 MinInGrpDefn[3]__2____31_____ <= 1
 MinInGrpDefn[3]__2____33_____ <= 1
 MinInGrpDefn[3]__2____35_____ <= 1
 MinInGrpDefn[3]__2____37_____ <= 1
 MinInGrpDefn[3]__2____39_____ <= 1
 MinInGrpDefn[3]__2____41_____ <= 1
 MinInGrpDefn[3]__2____43_____ <= 1
 MinInGrpDefn[3]__2____45_____ <= 1
 MinInGrpDefn[3]__2____47_____ <= 1
 MinInGrpDefn[3]__2____49_____ <= 1
 MinInGrpDefn[4]__2________ <= 1
 MinInGrpDefn[4]__2____5_____ <= 1
 MinInGrpDefn[4]__2____7_____ <= 1
 MinInGrpDefn[4]__2____9_____ <= 1
 MinInGrpDefn[4]__2____11_____ <= 1
 MinInGrpDefn[4]__2____13_____ <= 1
 MinInGrpDefn[4]__2____15_____ <= 1
 MinInGrpDefn[4]__2____17_____ <= 1
 MinInGrpDefn[4]__2____19_____ <= 1
 MinInGrpDefn[4]__2____21_____ <= 1
 MinInGrpDefn[4]__2____23_____ <= 1
 MinInGrpDefn[4]__2____25_____ <= 1
 MinInGrpDefn[4]__2____27_____ <= 1
 MinInGrpDefn[4]__2____31_____ <= 1
 MinInGrpDefn[4]__2____33_____ <= 1
 MinInGrpDefn[4]__2____35_____ <= 1
 MinInGrpDefn[4]__2____37_____ <= 1
 MinInGrpDefn[4]__2____39_____ <= 1
 MinInGrpDefn[4]__2____41_____ <= 1
 MinInGrpDefn[4]__2____43_____ <= 1
 MinInGrpDefn[4]__2____45_____ <= 1
 MinInGrpDefn[4]__2____47_____ <= 1
 MinInGrpDefn[4]__2____49_____ <= 1
 MinInGrpDefn[5]__2________ <= 1
 MinInGrpDefn[5]__2____5_____ <= 1
 MinInGrpDefn[5]__2____7_____ <= 1
 MinInGrpDefn[5]__2____9_____ <= 1
 MinInGrpDefn[5]__2____11_____ <= 1
 MinInGrpDefn[5]__2____13_____ <= 1
 MinInGrpDefn[5]__2____15_____ <= 1
 MinInGrpDefn[5]__2____17_____ <= 1
 MinInGrpDefn[5]__2____19_____ <= 1
 MinInGrpDefn[5]__2____21_____ <= 1
 MinInGrpDefn[5]__2____23_____ <= 1
 MinInGrpDefn[5]__2____25_____ <= 1
 MinInGrpDefn[5]__2____27_____ <= 1
 MinInGrpDefn[5]__2____31_____ <= 1
 MinInGrpDefn[5]__2____33_____ <= 1
 MinInGrpDefn[5]__2____35_____ <= 1
 MinInGrpDefn[5]__2____37_____ <= 1
 MinInGrpDefn[5]__2____39_____ <= 1
 MinInGrpDefn[5]__2____41_____ <= 1
 MinInGrpDefn[5]__2____43_____ <= 1
 MinInGrpDefn[5]__2____45_____ <= 1
 MinInGrpDefn[5]__2____47_____ <= 1
 MinInGrpDefn[5]__2____49_____ <= 1
 MinInGrpDefn[4]__2____50________ <= 1
 MinInGrpDefn[4]__2____50____2_____ <= 1
 MinInGrpDefn[3]__2____50________ <= 1
 MinInGrpDefn[3]__2____50____2_____ <= 1
 MinInGrpDefn[5]__2____50________ <= 1
 MinInGrpDefn[5]__2____50____2_____ <= 1
 MinInGrpDefn[2]__2____50________ <= 1
 MinInGrpDefn[2]__2____50____2_____ <= 1
 MinInGrpDefn[1]__2____50________ <= 1
 MinInGrpDefn[1]__2____50____2_____ <= 1
 MinInGrpDefn[4]__2____48________ <= 1
 MinInGrpDefn[4]__2____48____2_____ <= 1
 MinInGrpDefn[3]__2____48________ <= 1
 MinInGrpDefn[3]__2____48____2_____ <= 1
 MinInGrpDefn[5]__2____48________ <= 1
 MinInGrpDefn[5]__2____48____2_____ <= 1
 MinInGrpDefn[2]__2____48________ <= 1
 MinInGrpDefn[2]__2____48____2_____ <= 1
 MinInGrpDefn[1]__2____48________ <= 1
 MinInGrpDefn[1]__2____48____2_____ <= 1
 MinInGrpDefn[4]__2____46________ <= 1
 MinInGrpDefn[4]__2____46____2_____ <= 1
 MinInGrpDefn[3]__2____46________ <= 1
 MinInGrpDefn[3]__2____46____2_____ <= 1
 MinInGrpDefn[5]__2____46________ <= 1
 MinInGrpDefn[5]__2____46____2_____ <= 1
 MinInGrpDefn[2]__2____46________ <= 1
 MinInGrpDefn[2]__2____46____2_____ <= 1
 MinInGrpDefn[1]__2____46________ <= 1
 MinInGrpDefn[1]__2____46____2_____ <= 1
 MinInGrpDefn[4]__2____44________ <= 1
 MinInGrpDefn[4]__2____44____2_____ <= 1
 MinInGrpDefn[3]__2____44________ <= 1
 MinInGrpDefn[3]__2____44____2_____ <= 1
 MinInGrpDefn[5]__2____44________ <= 1
 MinInGrpDefn[5]__2____44____2_____ <= 1
 MinInGrpDefn[2]__2____44________ <= 1
 MinInGrpDefn[2]__2____44____2_____ <= 1
 MinInGrpDefn[1]__2____44________ <= 1
 MinInGrpDefn[1]__2____44____2_____ <= 1
 MinInGrpDefn[4]__2____42________ <= 1
 MinInGrpDefn[4]__2____42____2_____ <= 1
 MinInGrpDefn[3]__2____42________ <= 1
 MinInGrpDefn[3]__2____42____2_____ <= 1
 MinInGrpDefn[5]__2____42________ <= 1
 MinInGrpDefn[5]__2____42____2_____ <= 1
 MinInGrpDefn[2]__2____42________ <= 1
 MinInGrpDefn[2]__2____42____2_____ <= 1
 MinInGrpDefn[1]__2____42________ <= 1
 MinInGrpDefn[1]__2____42____2_____ <= 1
 MinInGrpDefn[4]__2____40________ <= 1
 MinInGrpDefn[4]__2____40____2_____ <= 1
 MinInGrpDefn[3]__2____40________ <= 1
 MinInGrpDefn[3]__2____40____2_____ <= 1
 MinInGrpDefn[5]__2____40________ <= 1
 MinInGrpDefn[5]__2____40____2_____ <= 1
 MinInGrpDefn[2]__2____40________ <= 1
 MinInGrpDefn[2]__2____40____2_____ <= 1
 MinInGrpDefn[1]__2____40________ <= 1
 MinInGrpDefn[1]__2____40____2_____ <= 1
 MinInGrpDefn[4]__2____38________ <= 1
 MinInGrpDefn[4]__2____38____2_____ <= 1
 MinInGrpDefn[3]__2____38________ <= 1
 MinInGrpDefn[3]__2____38____2_____ <= 1
 MinInGrpDefn[5]__2____38________ <= 1
 MinInGrpDefn[5]__2____38____2_____ <= 1
 MinInGrpDefn[2]__2____38________ <= 1
 MinInGrpDefn[2]__2____38____2_____ <= 1
 MinInGrpDefn[1]__2____38________ <= 1
 MinInGrpDefn[1]__2____38____2_____ <= 1
 MinInGrpDefn[4]__2____36________ <= 1
 MinInGrpDefn[4]__2____36____2_____ <= 1
 MinInGrpDefn[3]__2____36________ <= 1
 MinInGrpDefn[3]__2____36____2_____ <= 1
 MinInGrpDefn[5]__2____36________ <= 1
 MinInGrpDefn[5]__2____36____2_____ <= 1
 MinInGrpDefn[2]__2____36________ <= 1
 MinInGrpDefn[2]__2____36____2_____ <= 1
 MinInGrpDefn[1]__2____36________ <= 1
 MinInGrpDefn[1]__2____36____2_____ <= 1
 MinInGrpDefn[4]__2____34________ <= 1
 MinInGrpDefn[4]__2____34____2_____ <= 1
 MinInGrpDefn[3]__2____34________ <= 1
 MinInGrpDefn[3]__2____34____2_____ <= 1
 MinInGrpDefn[5]__2____34________ <= 1
 MinInGrpDefn[5]__2____34____2_____ <= 1
 MinInGrpDefn[2]__2____34________ <= 1
 MinInGrpDefn[2]__2____34____2_____ <= 1
 MinInGrpDefn[1]__2____34________ <= 1
 MinInGrpDefn[1]__2____34____2_____ <= 1
 MinInGrpDefn[4]__2____32________ <= 1
 MinInGrpDefn[4]__2____32____2_____ <= 1
 MinInGrpDefn[3]__2____32________ <= 1
 MinInGrpDefn[3]__2____32____2_____ <= 1
 MinInGrpDefn[5]__2____32________ <= 1
 MinInGrpDefn[5]__2____32____2_____ <= 1
 MinInGrpDefn[2]__2____32________ <= 1
 MinInGrpDefn[2]__2____32____2_____ <= 1
 MinInGrpDefn[1]__2____32________ <= 1
 MinInGrpDefn[1]__2____32____2_____ <= 1
 MinInGrpDefn[4]__2____28________ <= 1
 MinInGrpDefn[4]__2____28____2_____ <= 1
 MinInGrpDefn[3]__2____28________ <= 1
 MinInGrpDefn[3]__2____28____2_____ <= 1
 MinInGrpDefn[5]__2____28________ <= 1
 MinInGrpDefn[5]__2____28____2_____ <= 1
 MinInGrpDefn[2]__2____28________ <= 1
 MinInGrpDefn[2]__2____28____2_____ <= 1
 MinInGrpDefn[1]__2____28________ <= 1
 MinInGrpDefn[1]__2____28____2_____ <= 1
 MinInGrpDefn[4]__2____26________ <= 1
 MinInGrpDefn[4]__2____26____2_____ <= 1
 MinInGrpDefn[3]__2____26________ <= 1
 MinInGrpDefn[3]__2____26____2_____ <= 1
 MinInGrpDefn[5]__2____26________ <= 1
 MinInGrpDefn[5]__2____26____2_____ <= 1
 MinInGrpDefn[2]__2____26________ <= 1
 MinInGrpDefn[2]__2____26____2_____ <= 1
 MinInGrpDefn[1]__2____26________ <= 1
 MinInGrpDefn[1]__2____26____2_____ <= 1
 MinInGrpDefn[4]__2____24________ <= 1
 MinInGrpDefn[4]__2____24____2_____ <= 1
 MinInGrpDefn[3]__2____24________ <= 1
 MinInGrpDefn[3]__2____24____2_____ <= 1
 MinInGrpDefn[5]__2____24________ <= 1
 MinInGrpDefn[5]__2____24____2_____ <= 1
 MinInGrpDefn[2]__2____24________ <= 1
 MinInGrpDefn[2]__2____24____2_____ <= 1
 MinInGrpDefn[1]__2____24________ <= 1
 MinInGrpDefn[1]__2____24____2_____ <= 1
 MinInGrpDefn[4]__2____22________ <= 1
 MinInGrpDefn[4]__2____22____2_____ <= 1
 MinInGrpDefn[3]__2____22________ <= 1
 MinInGrpDefn[3]__2____22____2_____ <= 1
 MinInGrpDefn[5]__2____22________ <= 1
 MinInGrpDefn[5]__2____22____2_____ <= 1
 MinInGrpDefn[2]__2____22________ <= 1
 MinInGrpDefn[2]__2____22____2_____ <= 1
 MinInGrpDefn[1]__2____22________ <= 1
 MinInGrpDefn[1]__2____22____2_____ <= 1
 MinInGrpDefn[4]__2____20________ <= 1
 MinInGrpDefn[4]__2____20____2_____ <= 1
 MinInGrpDefn[3]__2____20________ <= 1
 MinInGrpDefn[3]__2____20____2_____ <= 1
 MinInGrpDefn[5]__2____20________ <= 1
 MinInGrpDefn[5]__2____20____2_____ <= 1
 MinInGrpDefn[2]__2____20________ <= 1
 MinInGrpDefn[2]__2____20____2_____ <= 1
 MinInGrpDefn[1]__2____20________ <= 1
 MinInGrpDefn[1]__2____20____2_____ <= 1
 MinInGrpDefn[4]__2____18________ <= 1
 MinInGrpDefn[4]__2____18____2_____ <= 1
 MinInGrpDefn[3]__2____18________ <= 1
 MinInGrpDefn[3]__2____18____2_____ <= 1
 MinInGrpDefn[5]__2____18________ <= 1
 MinInGrpDefn[5]__2____18____2_____ <= 1
 MinInGrpDefn[2]__2____18________ <= 1
 MinInGrpDefn[2]__2____18____2_____ <= 1
 MinInGrpDefn[1]__2____18________ <= 1
 MinInGrpDefn[1]__2____18____2_____ <= 1
 MinInGrpDefn[4]__2____16________ <= 1
 MinInGrpDefn[4]__2____16____2_____ <= 1
 MinInGrpDefn[3]__2____16________ <= 1
 MinInGrpDefn[3]__2____16____2_____ <= 1
 MinInGrpDefn[5]__2____16________ <= 1
 MinInGrpDefn[5]__2____16____2_____ <= 1
 MinInGrpDefn[2]__2____16________ <= 1
 MinInGrpDefn[2]__2____16____2_____ <= 1
 MinInGrpDefn[1]__2____16________ <= 1
 MinInGrpDefn[1]__2____16____2_____ <= 1
 MinInGrpDefn[4]__2____14________ <= 1
 MinInGrpDefn[4]__2____14____2_____ <= 1
 MinInGrpDefn[3]__2____14________ <= 1
 MinInGrpDefn[3]__2____14____2_____ <= 1
 MinInGrpDefn[5]__2____14________ <= 1
 MinInGrpDefn[5]__2____14____2_____ <= 1
 MinInGrpDefn[2]__2____14________ <= 1
 MinInGrpDefn[2]__2____14____2_____ <= 1
 MinInGrpDefn[1]__2____14________ <= 1
 MinInGrpDefn[1]__2____14____2_____ <= 1
 MinInGrpDefn[4]__2____12________ <= 1
 MinInGrpDefn[4]__2____12____2_____ <= 1
 MinInGrpDefn[3]__2____12________ <= 1
 MinInGrpDefn[3]__2____12____2_____ <= 1
 MinInGrpDefn[5]__2____12________ <= 1
 MinInGrpDefn[5]__2____12____2_____ <= 1
 MinInGrpDefn[2]__2____12________ <= 1
 MinInGrpDefn[2]__2____12____2_____ <= 1
 MinInGrpDefn[1]__2____12________ <= 1
 MinInGrpDefn[1]__2____12____2_____ <= 1
 MinInGrpDefn[4]__2____10________ <= 1
 MinInGrpDefn[4]__2____10____2_____ <= 1
 MinInGrpDefn[3]__2____10________ <= 1
 MinInGrpDefn[3]__2____10____2_____ <= 1
 MinInGrpDefn[5]__2____10________ <= 1
 MinInGrpDefn[5]__2____10____2_____ <= 1
 MinInGrpDefn[2]__2____10________ <= 1
 MinInGrpDefn[2]__2____10____2_____ <= 1
 MinInGrpDefn[1]__2____10________ <= 1
 MinInGrpDefn[1]__2____10____2_____ <= 1
 MinInGrpDefn[4]__2____6________ <= 1
 MinInGrpDefn[4]__2____6____2_____ <= 1
 MinInGrpDefn[3]__2____6________ <= 1
 MinInGrpDefn[3]__2____6____2_____ <= 1
 MinInGrpDefn[5]__2____6________ <= 1
 MinInGrpDefn[5]__2____6____2_____ <= 1
 MinInGrpDefn[2]__2____6________ <= 1
 MinInGrpDefn[2]__2____6____2_____ <= 1
 MinInGrpDefn[1]__2____6________ <= 1
 MinInGrpDefn[1]__2____6____2_____ <= 1
 MinInGrpDefn[4]__2____2________ <= 1
 MinInGrpDefn[4]__2____2____2_____ <= 1
 MinInGrpDefn[3]__2____2________ <= 1
 MinInGrpDefn[3]__2____2____2_____ <= 1
 MinInGrpDefn[5]__2____2________ <= 1
 MinInGrpDefn[5]__2____2____2_____ <= 1
 MinInGrpDefn[2]__2____2________ <= 1
 MinInGrpDefn[2]__2____2____2_____ <= 1
 MinInGrpDefn[1]__2____2________ <= 1
 MinInGrpDefn[1]__2____2____2_____ <= 1
 MinInGrpDefn[4]__2____8________ <= 1
 MinInGrpDefn[4]__2____8____2_____ <= 1
 MinInGrpDefn[3]__2____8________ <= 1
 MinInGrpDefn[3]__2____8____2_____ <= 1
 MinInGrpDefn[5]__2____8________ <= 1
 MinInGrpDefn[5]__2____8____2_____ <= 1
 MinInGrpDefn[2]__2____8________ <= 1
 MinInGrpDefn[2]__2____8____2_____ <= 1
 MinInGrpDefn[1]__2____8________ <= 1
 MinInGrpDefn[1]__2____8____2_____ <= 1
Generals
 Assign['SML']____ Assign['BCG']____ Assign['JAE']____ Assign['FTR']____
 Assign['JLN']____ Assign['SDN']____ Assign['BLE']____ Assign['ACI']____
 Assign['CWU']____ Assign['JWR']____ Assign['MBT']____ Assign['SCR']____
 Assign['RCR']____ Assign['RVY']____ Assign['KWN']____ Assign['PST']____
 Assign['VMA']____ Assign['RHA']____ Assign['FEE']____ Assign['RGE']____
 Assign['PSA']____ Assign['JMA']____ Assign['PSO']____ Assign['PSN']____
 Assign['JEN']____ Assign['RSN']____ MinInGrp____ MaxInGrp____
 MinType['dept','NNE']____ MinType['dept','NNW']____
 MinType['dept','WSW']____ MinType['dept','SSE']____
 MinType['loc','Springfield']____ MinType['loc','Urbana']____
 MinType['loc','Peoria']____ MinType['loc','Carbondale']____
 MinType['loc','Joliet']____ MinType['rate','A']____
 MinType['rate','B']____ MinType['title','Assistant']____
 MinType['title','Adjunct']____ MinType['title','Deputy']____
 MaxType['dept','NNE']____ MaxType['dept','NNW']____
 MaxType['dept','WSW']____ MaxType['dept','SSE']____
 MaxType['loc','Springfield']____ MaxType['loc','Urbana']____
 MaxType['loc','Peoria']____ MaxType['loc','Carbondale']____
 MaxType['loc','Joliet']____ MaxType['rate','A']____
 MaxType['rate','B']____ MaxType['title','Assistant']____
 MaxType['title','Adjunct']____ MaxType['title','Deputy']____
 MinInGrpDefn[1]______ MinInGrpDefn[2]______ MinInGrpDefn[3]______
 MinInGrpDefn[4]______ MinInGrpDefn[5]______
 MinTypeDefn[1,'dept','NNE']______ MinTypeDefn[1,'dept','NNW']______
 MinTypeDefn[1,'dept','WSW']______ MinTypeDefn[1,'dept','SSE']______
 MinTypeDefn[1,'loc','Springfield']______
 MinTypeDefn[1,'loc','Urbana']______ MinTypeDefn[1,'loc','Peoria']______
 MinTypeDefn[1,'loc','Carbondale']______
 MinTypeDefn[1,'loc','Joliet']______ MinTypeDefn[1,'rate','A']______
 MinTypeDefn[1,'rate','B']______ MinTypeDefn[1,'title','Assistant']______
 MinTypeDefn[1,'title','Adjunct']______
 MinTypeDefn[1,'title','Deputy']______ MinTypeDefn[2,'dept','NNE']______
 MinTypeDefn[2,'dept','NNW']______ MinTypeDefn[2,'dept','WSW']______
 MinTypeDefn[2,'dept','SSE']______ MinTypeDefn[2,'loc','Springfield']______
 MinTypeDefn[2,'loc','Urbana']______ MinTypeDefn[2,'loc','Peoria']______
 MinTypeDefn[2,'loc','Carbondale']______
 MinTypeDefn[2,'loc','Joliet']______ MinTypeDefn[2,'rate','A']______
 MinTypeDefn[2,'rate','B']______ MinTypeDefn[2,'title','Assistant']______
 MinTypeDefn[2,'title','Adjunct']______
 MinTypeDefn[2,'title','Deputy']______ MinTypeDefn[3,'dept','NNE']______
 MinTypeDefn[3,'dept','NNW']______ MinTypeDefn[3,'dept','WSW']______
 MinTypeDefn[3,'dept','SSE']______ MinTypeDefn[3,'loc','Springfield']______
 MinTypeDefn[3,'loc','Urbana']______ MinTypeDefn[3,'loc','Peoria']______
 MinTypeDefn[3,'loc','Carbondale']______
 MinTypeDefn[3,'loc','Joliet']______ MinTypeDefn[3,'rate','A']______
 MinTypeDefn[3,'rate','B']______ MinTypeDefn[3,'title','Assistant']______
 MinTypeDefn[3,'title','Adjunct']______
 MinTypeDefn[3,'title','Deputy']______ MinTypeDefn[4,'dept','NNE']______
 MinTypeDefn[4,'dept','NNW']______ MinTypeDefn[4,'dept','WSW']______
 MinTypeDefn[4,'dept','SSE']______ MinTypeDefn[4,'loc','Springfield']______
 MinTypeDefn[4,'loc','Urbana']______ MinTypeDefn[4,'loc','Peoria']______
 MinTypeDefn[4,'loc','Carbondale']______
 MinTypeDefn[4,'loc','Joliet']______ MinTypeDefn[4,'rate','A']______
 MinTypeDefn[4,'rate','B']______ MinTypeDefn[4,'title','Assistant']______
 MinTypeDefn[4,'title','Adjunct']______
 MinTypeDefn[4,'title','Deputy']______ MinTypeDefn[5,'dept','NNE']______
 MinTypeDefn[5,'dept','NNW']______ MinTypeDefn[5,'dept','WSW']______
 MinTypeDefn[5,'dept','SSE']______ MinTypeDefn[5,'loc','Springfield']______
 MinTypeDefn[5,'loc','Urbana']______ MinTypeDefn[5,'loc','Peoria']______
 MinTypeDefn[5,'loc','Carbondale']______
 MinTypeDefn[5,'loc','Joliet']______ MinTypeDefn[5,'rate','A']______
 MinTypeDefn[5,'rate','B']______ MinTypeDefn[5,'title','Assistant']______
 MinTypeDefn[5,'title','Adjunct']______
 MinTypeDefn[5,'title','Deputy']______ MinInGrpDefn[1]__2________
 MinInGrpDefn[1]__2____5_____ MinInGrpDefn[1]__2____7_____
 MinInGrpDefn[1]__2____9_____ MinInGrpDefn[1]__2____11_____
 MinInGrpDefn[1]__2____13_____ MinInGrpDefn[1]__2____15_____
 MinInGrpDefn[1]__2____17_____ MinInGrpDefn[1]__2____19_____
 MinInGrpDefn[1]__2____21_____ MinInGrpDefn[1]__2____23_____
 MinInGrpDefn[1]__2____25_____ MinInGrpDefn[1]__2____27_____
 MinInGrpDefn[1]__2____31_____ MinInGrpDefn[1]__2____33_____
 MinInGrpDefn[1]__2____35_____ MinInGrpDefn[1]__2____37_____
 MinInGrpDefn[1]__2____39_____ MinInGrpDefn[1]__2____41_____
 MinInGrpDefn[1]__2____43_____ MinInGrpDefn[1]__2____45_____
 MinInGrpDefn[1]__2____47_____ MinInGrpDefn[1]__2____49_____
 MinInGrpDefn[2]__2________ MinInGrpDefn[2]__2____5_____
 MinInGrpDefn[2]__2____7_____ MinInGrpDefn[2]__2____9_____
 MinInGrpDefn[2]__2____11_____ MinInGrpDefn[2]__2____13_____
 MinInGrpDefn[2]__2____15_____ MinInGrpDefn[2]__2____17_____
 MinInGrpDefn[2]__2____19_____ MinInGrpDefn[2]__2____21_____
 MinInGrpDefn[2]__2____23_____ MinInGrpDefn[2]__2____25_____
 MinInGrpDefn[2]__2____27_____ MinInGrpDefn[2]__2____31_____
 MinInGrpDefn[2]__2____33_____ MinInGrpDefn[2]__2____35_____
 MinInGrpDefn[2]__2____37_____ MinInGrpDefn[2]__2____39_____
 MinInGrpDefn[2]__2____41_____ MinInGrpDefn[2]__2____43_____
 MinInGrpDefn[2]__2____45_____ MinInGrpDefn[2]__2____47_____
 MinInGrpDefn[2]__2____49_____ MinInGrpDefn[3]__2________
 MinInGrpDefn[3]__2____5_____ MinInGrpDefn[3]__2____7_____
 MinInGrpDefn[3]__2____9_____ MinInGrpDefn[3]__2____11_____
 MinInGrpDefn[3]__2____13_____ MinInGrpDefn[3]__2____15_____
 MinInGrpDefn[3]__2____17_____ MinInGrpDefn[3]__2____19_____
 MinInGrpDefn[3]__2____21_____ MinInGrpDefn[3]__2____23_____
 MinInGrpDefn[3]__2____25_____ MinInGrpDefn[3]__2____27_____
 MinInGrpDefn[3]__2____31_____ MinInGrpDefn[3]__2____33_____
 MinInGrpDefn[3]__2____35_____ MinInGrpDefn[3]__2____37_____
 MinInGrpDefn[3]__2____39_____ MinInGrpDefn[3]__2____41_____
 MinInGrpDefn[3]__2____43_____ MinInGrpDefn[3]__2____45_____
 MinInGrpDefn[3]__2____47_____ MinInGrpDefn[3]__2____49_____
 MinInGrpDefn[4]__2________ MinInGrpDefn[4]__2____5_____
 MinInGrpDefn[4]__2____7_____ MinInGrpDefn[4]__2____9_____
 MinInGrpDefn[4]__2____11_____ MinInGrpDefn[4]__2____13_____
 MinInGrpDefn[4]__2____15_____ MinInGrpDefn[4]__2____17_____
 MinInGrpDefn[4]__2____19_____ MinInGrpDefn[4]__2____21_____
 MinInGrpDefn[4]__2____23_____ MinInGrpDefn[4]__2____25_____
 MinInGrpDefn[4]__2____27_____ MinInGrpDefn[4]__2____31_____
 MinInGrpDefn[4]__2____33_____ MinInGrpDefn[4]__2____35_____
 MinInGrpDefn[4]__2____37_____ MinInGrpDefn[4]__2____39_____
 MinInGrpDefn[4]__2____41_____ MinInGrpDefn[4]__2____43_____
 MinInGrpDefn[4]__2____45_____ MinInGrpDefn[4]__2____47_____
 MinInGrpDefn[4]__2____49_____ MinInGrpDefn[5]__2________
 MinInGrpDefn[5]__2____5_____ MinInGrpDefn[5]__2____7_____
 MinInGrpDefn[5]__2____9_____ MinInGrpDefn[5]__2____11_____
 MinInGrpDefn[5]__2____13_____ MinInGrpDefn[5]__2____15_____
 MinInGrpDefn[5]__2____17_____ MinInGrpDefn[5]__2____19_____
 MinInGrpDefn[5]__2____21_____ MinInGrpDefn[5]__2____23_____
 MinInGrpDefn[5]__2____25_____ MinInGrpDefn[5]__2____27_____
 MinInGrpDefn[5]__2____31_____ MinInGrpDefn[5]__2____33_____
 MinInGrpDefn[5]__2____35_____ MinInGrpDefn[5]__2____37_____
 MinInGrpDefn[5]__2____39_____ MinInGrpDefn[5]__2____41_____
 MinInGrpDefn[5]__2____43_____ MinInGrpDefn[5]__2____45_____
 MinInGrpDefn[5]__2____47_____ MinInGrpDefn[5]__2____49_____
 MinInGrpDefn[4]__2____50________ MinInGrpDefn[4]__2____50____2_____
 MinInGrpDefn[3]__2____50________ MinInGrpDefn[3]__2____50____2_____
 MinInGrpDefn[5]__2____50________ MinInGrpDefn[5]__2____50____2_____
 MinInGrpDefn[2]__2____50________ MinInGrpDefn[2]__2____50____2_____
 MinInGrpDefn[1]__2____50________ MinInGrpDefn[1]__2____50____2_____
 MinInGrpDefn[4]__2____48________ MinInGrpDefn[4]__2____48____2_____
 MinInGrpDefn[3]__2____48________ MinInGrpDefn[3]__2____48____2_____
 MinInGrpDefn[5]__2____48________ MinInGrpDefn[5]__2____48____2_____
 MinInGrpDefn[2]__2____48________ MinInGrpDefn[2]__2____48____2_____
 MinInGrpDefn[1]__2____48________ MinInGrpDefn[1]__2____48____2_____
 MinInGrpDefn[4]__2____46________ MinInGrpDefn[4]__2____46____2_____
 MinInGrpDefn[3]__2____46________ MinInGrpDefn[3]__2____46____2_____
 MinInGrpDefn[5]__2____46________ MinInGrpDefn[5]__2____46____2_____
 MinInGrpDefn[2]__2____46________ MinInGrpDefn[2]__2____46____2_____
 MinInGrpDefn[1]__2____46________ MinInGrpDefn[1]__2____46____2_____
 MinInGrpDefn[4]__2____44________ MinInGrpDefn[4]__2____44____2_____
 MinInGrpDefn[3]__2____44________ MinInGrpDefn[3]__2____44____2_____
 MinInGrpDefn[5]__2____44________ MinInGrpDefn[5]__2____44____2_____
 MinInGrpDefn[2]__2____44________ MinInGrpDefn[2]__2____44____2_____
 MinInGrpDefn[1]__2____44________ MinInGrpDefn[1]__2____44____2_____
 MinInGrpDefn[4]__2____42________ MinInGrpDefn[4]__2____42____2_____
 MinInGrpDefn[3]__2____42________ MinInGrpDefn[3]__2____42____2_____
 MinInGrpDefn[5]__2____42________ MinInGrpDefn[5]__2____42____2_____
 MinInGrpDefn[2]__2____42________ MinInGrpDefn[2]__2____42____2_____
 MinInGrpDefn[1]__2____42________ MinInGrpDefn[1]__2____42____2_____
 MinInGrpDefn[4]__2____40________ MinInGrpDefn[4]__2____40____2_____
 MinInGrpDefn[3]__2____40________ MinInGrpDefn[3]__2____40____2_____
 MinInGrpDefn[5]__2____40________ MinInGrpDefn[5]__2____40____2_____
 MinInGrpDefn[2]__2____40________ MinInGrpDefn[2]__2____40____2_____
 MinInGrpDefn[1]__2____40________ MinInGrpDefn[1]__2____40____2_____
 MinInGrpDefn[4]__2____38________ MinInGrpDefn[4]__2____38____2_____
 MinInGrpDefn[3]__2____38________ MinInGrpDefn[3]__2____38____2_____
 MinInGrpDefn[5]__2____38________ MinInGrpDefn[5]__2____38____2_____
 MinInGrpDefn[2]__2____38________ MinInGrpDefn[2]__2____38____2_____
 MinInGrpDefn[1]__2____38________ MinInGrpDefn[1]__2____38____2_____
 MinInGrpDefn[4]__2____36________ MinInGrpDefn[4]__2____36____2_____
 MinInGrpDefn[3]__2____36________ MinInGrpDefn[3]__2____36____2_____
 MinInGrpDefn[5]__2____36________ MinInGrpDefn[5]__2____36____2_____
 MinInGrpDefn[2]__2____36________ MinInGrpDefn[2]__2____36____2_____
 MinInGrpDefn[1]__2____36________ MinInGrpDefn[1]__2____36____2_____
 MinInGrpDefn[4]__2____34________ MinInGrpDefn[4]__2____34____2_____
 MinInGrpDefn[3]__2____34________ MinInGrpDefn[3]__2____34____2_____
 MinInGrpDefn[5]__2____34________ MinInGrpDefn[5]__2____34____2_____
 MinInGrpDefn[2]__2____34________ MinInGrpDefn[2]__2____34____2_____
 MinInGrpDefn[1]__2____34________ MinInGrpDefn[1]__2____34____2_____
 MinInGrpDefn[4]__2____32________ MinInGrpDefn[4]__2____32____2_____
 MinInGrpDefn[3]__2____32________ MinInGrpDefn[3]__2____32____2_____
 MinInGrpDefn[5]__2____32________ MinInGrpDefn[5]__2____32____2_____
 MinInGrpDefn[2]__2____32________ MinInGrpDefn[2]__2____32____2_____
 MinInGrpDefn[1]__2____32________ MinInGrpDefn[1]__2____32____2_____
 MinInGrpDefn[4]__2____28________ MinInGrpDefn[4]__2____28____2_____
 MinInGrpDefn[3]__2____28________ MinInGrpDefn[3]__2____28____2_____
 MinInGrpDefn[5]__2____28________ MinInGrpDefn[5]__2____28____2_____
 MinInGrpDefn[2]__2____28________ MinInGrpDefn[2]__2____28____2_____
 MinInGrpDefn[1]__2____28________ MinInGrpDefn[1]__2____28____2_____
 MinInGrpDefn[4]__2____26________ MinInGrpDefn[4]__2____26____2_____
 MinInGrpDefn[3]__2____26________ MinInGrpDefn[3]__2____26____2_____
 MinInGrpDefn[5]__2____26________ MinInGrpDefn[5]__2____26____2_____
 MinInGrpDefn[2]__2____26________ MinInGrpDefn[2]__2____26____2_____
 MinInGrpDefn[1]__2____26________ MinInGrpDefn[1]__2____26____2_____
 MinInGrpDefn[4]__2____24________ MinInGrpDefn[4]__2____24____2_____
 MinInGrpDefn[3]__2____24________ MinInGrpDefn[3]__2____24____2_____
 MinInGrpDefn[5]__2____24________ MinInGrpDefn[5]__2____24____2_____
 MinInGrpDefn[2]__2____24________ MinInGrpDefn[2]__2____24____2_____
 MinInGrpDefn[1]__2____24________ MinInGrpDefn[1]__2____24____2_____
 MinInGrpDefn[4]__2____22________ MinInGrpDefn[4]__2____22____2_____
 MinInGrpDefn[3]__2____22________ MinInGrpDefn[3]__2____22____2_____
 MinInGrpDefn[5]__2____22________ MinInGrpDefn[5]__2____22____2_____
 MinInGrpDefn[2]__2____22________ MinInGrpDefn[2]__2____22____2_____
 MinInGrpDefn[1]__2____22________ MinInGrpDefn[1]__2____22____2_____
 MinInGrpDefn[4]__2____20________ MinInGrpDefn[4]__2____20____2_____
 MinInGrpDefn[3]__2____20________ MinInGrpDefn[3]__2____20____2_____
 MinInGrpDefn[5]__2____20________ MinInGrpDefn[5]__2____20____2_____
 MinInGrpDefn[2]__2____20________ MinInGrpDefn[2]__2____20____2_____
 MinInGrpDefn[1]__2____20________ MinInGrpDefn[1]__2____20____2_____
 MinInGrpDefn[4]__2____18________ MinInGrpDefn[4]__2____18____2_____
 MinInGrpDefn[3]__2____18________ MinInGrpDefn[3]__2____18____2_____
 MinInGrpDefn[5]__2____18________ MinInGrpDefn[5]__2____18____2_____
 MinInGrpDefn[2]__2____18________ MinInGrpDefn[2]__2____18____2_____
 MinInGrpDefn[1]__2____18________ MinInGrpDefn[1]__2____18____2_____
 MinInGrpDefn[4]__2____16________ MinInGrpDefn[4]__2____16____2_____
 MinInGrpDefn[3]__2____16________ MinInGrpDefn[3]__2____16____2_____
 MinInGrpDefn[5]__2____16________ MinInGrpDefn[5]__2____16____2_____
 MinInGrpDefn[2]__2____16________ MinInGrpDefn[2]__2____16____2_____
 MinInGrpDefn[1]__2____16________ MinInGrpDefn[1]__2____16____2_____
 MinInGrpDefn[4]__2____14________ MinInGrpDefn[4]__2____14____2_____
 MinInGrpDefn[3]__2____14________ MinInGrpDefn[3]__2____14____2_____
 MinInGrpDefn[5]__2____14________ MinInGrpDefn[5]__2____14____2_____
 MinInGrpDefn[2]__2____14________ MinInGrpDefn[2]__2____14____2_____
 MinInGrpDefn[1]__2____14________ MinInGrpDefn[1]__2____14____2_____
 MinInGrpDefn[4]__2____12________ MinInGrpDefn[4]__2____12____2_____
 MinInGrpDefn[3]__2____12________ MinInGrpDefn[3]__2____12____2_____
 MinInGrpDefn[5]__2____12________ MinInGrpDefn[5]__2____12____2_____
 MinInGrpDefn[2]__2____12________ MinInGrpDefn[2]__2____12____2_____
 MinInGrpDefn[1]__2____12________ MinInGrpDefn[1]__2____12____2_____
 MinInGrpDefn[4]__2____10________ MinInGrpDefn[4]__2____10____2_____
 MinInGrpDefn[3]__2____10________ MinInGrpDefn[3]__2____10____2_____
 MinInGrpDefn[5]__2____10________ MinInGrpDefn[5]__2____10____2_____
 MinInGrpDefn[2]__2____10________ MinInGrpDefn[2]__2____10____2_____
 MinInGrpDefn[1]__2____10________ MinInGrpDefn[1]__2____10____2_____
 MinInGrpDefn[4]__2____6________ MinInGrpDefn[4]__2____6____2_____
 MinInGrpDefn[3]__2____6________ MinInGrpDefn[3]__2____6____2_____
 MinInGrpDefn[5]__2____6________ MinInGrpDefn[5]__2____6____2_____
 MinInGrpDefn[2]__2____6________ MinInGrpDefn[2]__2____6____2_____
 MinInGrpDefn[1]__2____6________ MinInGrpDefn[1]__2____6____2_____
 MinInGrpDefn[4]__2____2________ MinInGrpDefn[4]__2____2____2_____
 MinInGrpDefn[3]__2____2________ MinInGrpDefn[3]__2____2____2_____
 MinInGrpDefn[5]__2____2________ MinInGrpDefn[5]__2____2____2_____
 MinInGrpDefn[2]__2____2________ MinInGrpDefn[2]__2____2____2_____
 MinInGrpDefn[1]__2____2________ MinInGrpDefn[1]__2____2____2_____
 MinInGrpDefn[4]__2____8________ MinInGrpDefn[4]__2____8____2_____
 MinInGrpDefn[3]__2____8________ MinInGrpDefn[3]__2____8____2_____
 MinInGrpDefn[5]__2____8________ MinInGrpDefn[5]__2____8____2_____
 MinInGrpDefn[2]__2____8________ MinInGrpDefn[2]__2____8____2_____
 MinInGrpDefn[1]__2____8________ MinInGrpDefn[1]__2____8____2_____
End
