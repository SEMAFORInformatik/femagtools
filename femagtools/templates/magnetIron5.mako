
m.magn_rad       = da2/2
m.yoke_rad       = dy2/2


	 CALL L_GETFLOAT("m.magn_di_ra",magn_di_ra,IRET)
	 CALL L_GETFLOAT("m.gap_ma_iron",gap_ma_iron,IRET)
	 CALL L_GETFLOAT("m.iron_bfe",iron_bfe,IRET)
	 CALL L_GETFLOAT("m.air_space_h",air_space_h,IRET)
	 CALL L_GETFLOAT("m.corner_r",corner_r,IRET)
	 CALL L_GETFLOAT("m.air_sp_ori",air_sp_ori,IRET)
	 CALL L_GETFLOAT("m.magn_num",magn_num,IRET)
	 CALL L_GETFLOAT("m.iron_shape",iron_shape,IRET)
	 CALL L_GETFLOAT("m.air_space_b",air_space_b,IRET)
	 CALL L_GETFLOAT("m.magn_width",magn_width,IRET)
	 CALL L_GETFLOAT("m.nodedist",nodedist,IRET)
	 CALL L_GETFLOAT("m.airgap",airgap,IRET)

m.magn_height     =    ${model.get(['magnet','magnetIron5', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetIron5', 'magn_width'])*1e3}
m.gap_ma_iron     =    ${model.get(['magnet','magnetIron5', 'gap_ma_iron'])*1e3}
m.iron_bfe        =    ${model.get(['magnet','magnetIron5', 'iron_bfe'])*1e3}
m.air_space_h     =    ${model.get(['magnet','magnetIron5', 'air_space_h'])*1e3}
m.corner_r        =    ${model.get(['magnet','magnetIron5', 'corner_r'])*1e3}
m.air_sp_ori      =    ${model.get(['magnet','magnetIron5', 'air_sp_ori'])}
m.magn_num        =    ${model.get(['magnet','magnetIron5', 'magn_num'])}
m.iron_shape      =    ${model.get(['magnet','magnetIron5', 'iron_shape'])*1e3}
m.air_space_b     =    ${model.get(['magnet','magnetIron5', 'air_space_b'])*1e3}
m.magn_di_ra      =    ${model.get(['magnet','magnetIron5', 'magn_di_ra'])*1e3}

m.zeroangl        =     0.0

m.mcvkey_yoke     =   '${model.magnet.get('mcvkey_yoke','dummy')}'
m.mcvkey_mshaft   =   '${model.magnet.get('mcvkey_mshaft','dummy')}'
m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet Iron 5")
