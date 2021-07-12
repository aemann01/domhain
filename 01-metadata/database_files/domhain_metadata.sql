-- -----------------------------------------------------
-- Schema domhain
-- -----------------------------------------------------
SET GLOBAL local_infile=1;
SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';

-- -----------------------------------------------------
-- remove database if it already exists
-- -----------------------------------------------------
DROP DATABASE IF EXISTS domhain;

-- -----------------------------------------------------
-- create new database
-- -----------------------------------------------------
CREATE DATABASE domhain;
USE domhain;

-- -----------------------------------------------------
-- Table domain.Demographic
-- -----------------------------------------------------
CREATE TABLE Demographic (
	study_id VARCHAR(255) NOT NULL,
	dob DATE,
	visit_num VARCHAR(255),
	age_y INT,
	age_m INT,
	sex VARCHAR(255),
	study_group VARCHAR(255),
	PRIMARY KEY (study_id)
);

-- -----------------------------------------------------
-- Table domain.Sample_manifest
-- -----------------------------------------------------
CREATE TABLE Sample_manifest (
	manifest_id VARCHAR(255) NOT NULL,
	study_id VARCHAR(255),
	visit_num INT,
	facility VARCHAR(255),
	sample_type VARCHAR(255),
	tube_type VARCHAR(255),
	collection_date DATE,
	collection_time TIME,
	aliquot_number INT,
	aliquot_type VARCHAR(255),
	aliquot_code VARCHAR(255),
	initial_ml VARCHAR(255),
	current_ml VARCHAR(255),
	box_num INT,
	row_num VARCHAR(255),
	col_num INT,
	hiv_status VARCHAR(255),
	PRIMARY KEY (manifest_id),
	FOREIGN KEY (study_id)
		REFERENCES Demographic(study_id)
);

-- -----------------------------------------------------
-- Table domain.Breast_feeding
-- -----------------------------------------------------
CREATE TABLE Breast_feeding (
	study_id VARCHAR(255) NOT NULL,
	breast_fed VARCHAR(255),
	breast_min INT,
	breast_hrs INT,
	breast_days INT,
	colostrum INT,
	water_before_breast_fed INT,
	infant_formula INT,
	drink_fluids INT,
	stop_breast_fed INT,
	why_not_breast_fed INT,
	PRIMARY KEY (study_id),
	FOREIGN KEY (study_id)
		REFERENCES Demographic(study_id)
);

-- -----------------------------------------------------
-- Table domain.Diet
-- -----------------------------------------------------
CREATE TABLE Diet (
	study_id VARCHAR(255) NOT NULL,
	visit_id INT,
	food VARCHAR(255) NOT NULL,
	is_food_taken INT,
	how_often INT,
	last_24hrs_frequency INT,
	PRIMARY KEY (study_id, food), 
	FOREIGN KEY (study_id)
		REFERENCES Demographic(study_id)
);

-- -----------------------------------------------------
-- Table domain.Diet_food_taken
-- -----------------------------------------------------
CREATE TABLE Diet_food_taken (
	study_id VARCHAR(255) NOT NULL,
	visit_num INT,
	antibiotics_syrup INT,
	beans INT,
	biscuits_cookies INT,
	bread INT,
	cake_buns_puffpuff INT,
	eba_amala_fufu INT,
	eggs INT,
	fizzy_drinks INT,
	fura INT,
	ice_cream INT,
	indomie_noodles INT,
	juices INT,
	maize_other_corn_meal INT,
	meat_fish_chicken INT,
	multivitamins_syrup INT,
	other_solids INT,
	pap_ogi_akamu INT,
	rice INT,
	soya_milk INT,
	sugar_or_glucose_water INT,
	sugar_salt_soin_ors INT,
	sugary_liquids INT,
	sweets_chocolate INT,
	tea_chocolate_drink INT,
	tinned_powdered_milk INT,
	vegetables_fruits INT,
	yam_yam_pottage INT,
	PRIMARY KEY (study_id),
	FOREIGN KEY (study_id)
		REFERENCES Demographic(study_id)
);

-- -----------------------------------------------------
-- Table domain.Diet_24hrs
-- -----------------------------------------------------
CREATE TABLE Diet_24hrs (
	study_id VARCHAR(255) NOT NULL,
	visit_num INT,
	antibiotics_syrup INT,
	beans INT,
	biscuits_cookies INT,
	bread INT,
	cake_buns_puffpuff INT,
	eba_amala_fufu INT,
	eggs INT,
	fizzy_drinks INT,
	fura INT,
	ice_cream INT,
	indomie_noodles INT,
	juices INT,
	maize_other_corn_meal INT,
	meat_fish_chicken INT,
	multivitamins_syrup INT,
	other_solids INT,
	pap_ogi_akamu INT,
	rice INT,
	soya_milk INT,
	sugar_or_glucose_water INT,
	sugar_salt_soin_ors INT,
	sugary_liquids INT,
	sweets_chocolate INT,
	tea_chocolate_drink INT,
	tinned_powdered_milk INT,
	vegetables_fruits INT,
	yam_yam_pottage INT,
	PRIMARY KEY (study_id),
	FOREIGN KEY (study_id)
		REFERENCES Demographic(study_id)
);

-- -----------------------------------------------------
-- Table domain.Clinical
-- -----------------------------------------------------
CREATE TABLE Clinical (
	study_id VARCHAR(255) NOT NULL,
	assign_date DATE,
	hiv_diagnosis DATE,
	delivery_method INT,
	membrane_rupture INT,
	gestational_age INT,
	birth_weight INT,
	pre_mature INT,
	arv_at_birth INT,
	art_during_pregnancy INT,
	definitive_diagnosis INT,
	last_pcr_or_antibody_test DATE,
	pcr_or_antibody_result INT,
	viral_load INT,
	aids_symptoms INT,
	symptom_type VARCHAR(255),
	pneumonia INT,
	no_medication INT,
	sickle_cell INT,
	measles INT,
	upper_respiratory_tract INT,
	other_medication VARCHAR(255),
	cd4_percent INT,
	cd4_per_mm3 INT,
	current_viral_load INT,
	current_weight INT,
	current_height INT,
	first_regimen_arv VARCHAR(255),
	first_dosage_form DATE,
	second_regimen_arv VARCHAR(255),
	second_dosage_form DATE,
	fixed_dose VARCHAR(255),
	fixed_dose_form VARCHAR(255),
	fixed_dose_am INT,
	fixed_dose_pm INT,
	fixed_treatment_duration VARCHAR(255),
	single_drugs VARCHAR(255),
	single_drugs_form VARCHAR(255),
	single_drugs_am INT,
	single_drugs_pm INT,
	single_drugs_duration VARCHAR(255),
	oi_drugs VARCHAR(255),
	oi_drugs_form VARCHAR(255),
	oi_drugs_am INT,
	oi_drugs_pm INT,
	oi_drugs_duration VARCHAR(255),
	staff_name VARCHAR(255),
	data_tech VARCHAR(255),
	PRIMARY KEY (study_id), 
	FOREIGN KEY (study_id)
		REFERENCES Demographic(study_id)
);

-- -----------------------------------------------------
-- Table domain.Lab
-- -----------------------------------------------------
CREATE TABLE Lab (
	sample_id VARCHAR(255) NOT NULL,
	study_id VARCHAR(255) NOT NULL,
	plaque_amount INT,
	ul_for_extraction INT,
	extraction_date DATE,
	dna_concentration_ngul FLOAT(3,2),
	PRIMARY KEY (sample_id),
	FOREIGN KEY (study_id)
		REFERENCES Demographic(study_id)
);

-- -----------------------------------------------------
-- Load data files
-- -----------------------------------------------------
-- NOTE: if you have a windows computer you'll need to change the carriage return (\r) to newline (\n)
LOAD DATA LOCAL INFILE 'sample_manifest.csv' INTO TABLE Sample_manifest FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'breast_feeding.csv' INTO TABLE Breast_feeding FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'demographic.csv' INTO TABLE Demographic FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'diet.csv' INTO TABLE Diet FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'diet_24hrs.csv' INTO TABLE Diet_24hrs FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'diet_food_taken.csv' INTO TABLE Diet_food_taken FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'clinical.csv' INTO TABLE Clinical FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'lab.csv' INTO TABLE Lab FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r' IGNORE 1 ROWS;