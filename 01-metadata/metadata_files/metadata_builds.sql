-- -----------------------------------------------------
-- DOMHAIN METADATA QUERYS
-- -----------------------------------------------------

-- Table containing dietary intake of sweets, medicine, or fermented foods in the last 24hrs before collection
CREATE TABLE Sweets_fermented_foods 
	SELECT My_samples.manifest_id, 
	My_samples.visit_num, 
	Demographic.study_id, 
	My_samples.aliquot_type,
	Demographic.sex, 
	Demographic.age_y, 
	Demographic.study_group,
	Diet_24hrs.fizzy_drinks,
	Diet_24hrs.antibiotics_syrup,
	Diet_24hrs.biscuits_cookies,
	Diet_24hrs.cake_buns_puffpuff,
	Diet_24hrs.sugar_or_glucose_water,
	Diet_24hrs.sugar_salt_soin_ors,
	Diet_24hrs.sugary_liquids,
	Diet_24hrs.sweets_chocolate,
	Diet_24hrs.fura,
	Diet_24hrs.eba_amala_fufu,
	Diet_24hrs.pap_ogi_akamu,
	Diet_24hrs.tea_chocolate_drink 
	FROM My_samples,Demographic,Diet_24hrs 
	WHERE My_samples.study_id=Demographic.study_id AND My_samples.study_id=Diet_24hrs.study_id;
-- save to file
SELECT * 
	FROM Sweets_fermented_foods 
	ORDER BY study_id INTO OUTFILE '/Users/mann/github/domhain/01-metadata/metadata_files/sweets_fermented_foods.txt' 
	FIELDS TERMINATED BY '\t' 
	LINES TERMINATED BY '\n';
