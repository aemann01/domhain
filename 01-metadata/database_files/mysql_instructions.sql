-- -----------------------------------------------------
-- MYSQL QUICK AND DIRTY TUTORIAL
-- -----------------------------------------------------

-- -----------------------------------------------------
-- SOME NOTES BEFORE YOU START
-- -----------------------------------------------------
-- MySQL doesn't like the , or & in the domhain sample names, be sure to change them to '-' in all of your self-made metadata files before starting
-- This tutorial is based on Mac OS with homebrew installed, setting MySQL up on a different OS will be different

-- -----------------------------------------------------
-- INSTALLATION
-- -----------------------------------------------------
-- install with homebrew
brew install mysql 
-- once mysql is installed you'll have to set up a secure installation
-- this is where you will be prompted to install a root password for your local mysql server
mysql_secure_installation
-- by default mysql won't let you easily read or write files outside of 'secure' folders, to disable this:
-- NOTE: be sure to check that this file /usr/local/etc/my.cnf exists before executing the command below
echo "[mysqld]\nsecure_file_priv\t\t= ''\n" | sudo tee /usr/local/etc/my.cnf

-- -----------------------------------------------------
-- LOAD DATABASE AND START LOCAL SERVER
-- -----------------------------------------------------
-- load database -- only will need to do this once unless updating part of the database
mysql -u root -p < domhain_metadata.sql
-- login to localhost 
mysql -u root -p

-- -----------------------------------------------------
-- GENERAL NAVIGATION
-- -----------------------------------------------------
-- what databases do we have available?
SHOW DATABASES;
-- besides the domhain database there will be a few others that come with your local install
-- load the domhain database
USE domhain;
-- within the domhain database, what tables are available?
SHOW TABLES;

-- -----------------------------------------------------
-- NAVIGATING TABLES
-- -----------------------------------------------------
-- see all contents of a table
SELECT * FROM Diet;
-- summary of table structure
DESCRIBE Demographic;
-- OR
SHOW COLUMNS 
	FROM Lab;
-- select a single column from table
SELECT dob
	FROM Demographic;
-- show first 10 lines of table
SELECT * FROM Sample_manifest 
LIMIT 10;
-- select all rows from demographic table where sex = female
SELECT * 
	FROM Demographic WHERE sex="F";
-- two conditions
SELECT * 
	FROM Demographic 
	WHERE sex="F" 
	AND study_group="HUU";
-- three condtions
SELECT * 
	FROM Demographic 
	WHERE sex="F" 
	AND study_group="HUU" 
	AND age_y >= "10";
-- or statement
SELECT * 
	FROM Demographic 
	WHERE sex="F" 
	AND (study_group="HUU" OR study_group="HEU");

-- -----------------------------------------------------
-- JOINING AND MANIPULATING TABLES
-- -----------------------------------------------------
-- this is where the real power of mysql lies, because tables are relational, you can easily query one from another and build custom metadata files
-- join two tables with shared unique values in columns
SELECT * 
	FROM Sample_manifest 
	JOIN Lab 
	ON Sample_manifest.manifest_id = Lab.sample_id;
-- another example
SELECT * 
	FROM Breast_feeding 
	JOIN Demographic 
	ON Breast_feeding.study_id = Demographic.study_id;

-- -----------------------------------------------------
-- CREATING YOUR OWN TABLES
-- -----------------------------------------------------
-- ok so now say you want to subset your metadata based on your sample sheet
-- remove a table (the first time you run this you'll get a 'this table doesn't exist' message - no worries)
DROP TABLE Sample_sheet;
-- first need to make a blank table to populate (can't use "index" as column name)
CREATE TABLE Sample_sheet (
	sample_id VARCHAR(255) NOT NULL,
	sample_name VARCHAR(255),
	sample_plate INT,
	sample_well INT,
	i5_index_id VARCHAR(255),
	index1 VARCHAR(255),
	i7_index_id VARCHAR(255),
	index2 VARCHAR(255),
	sample_project VARCHAR(255),
	description VARCHAR(255),
	PRIMARY KEY (sample_id) 
);
-- now load your csv file into the newly created table
LOAD DATA LOCAL INFILE '/Users/mann/github/domhain/01-metadata/SampleSheet.csv'
	INTO TABLE Sample_sheet
	FIELDS TERMINATED BY ','
	LINES TERMINATED BY '\n'
	IGNORE 18 LINES
;
-- not loading properly? check warnings
SHOW WARNINGS;
-- check your file
SELECT * FROM Sample_sheet LIMIT 5;
-- subsample the sample manifest 
SELECT * 
	FROM Sample_manifest 
	JOIN Sample_sheet 
	ON Sample_manifest.manifest_id = Sample_sheet.sample_id;
-- create a new table from the join
CREATE TABLE My_samples
SELECT * 
	FROM Sample_manifest 
	JOIN Sample_sheet 
	ON Sample_manifest.manifest_id = Sample_sheet.sample_id;
-- see head of new table
SELECT * FROM My_samples LIMIT 5;
-- how many unique samples IDs are in our data?
SELECT COUNT(DISTINCT(manifest_id)) FROM My_samples;
-- how many individuals does this represent?
SELECT COUNT(DISTINCT(study_id)) FROM My_samples;
-- how many indiviuals fall under each hiv status category?
SELECT hiv_status, COUNT(hiv_status) FROM My_samples GROUP BY hiv_status;
-- you can now use your data to get other information from other relational databases
-- simple inner join of your samples and the demographic table
SELECT * 
	FROM My_samples,Demographic 
	WHERE My_samples.study_id=Demographic.study_id;
-- ok so now let's say we want to clean this up and output it to a file
-- first let's see which columns we want to keep in the merged table
DESCRIBE Demographic;
-- +-------------+--------------+------+-----+---------+-------+
-- | Field       | Type         | Null | Key | Default | Extra |
-- +-------------+--------------+------+-----+---------+-------+
-- | study_id    | varchar(255) | NO   | PRI | NULL    |       |
-- | dob         | date         | YES  |     | NULL    |       |
-- | visit_num   | varchar(255) | YES  |     | NULL    |       |
-- | age_y       | int          | YES  |     | NULL    |       |
-- | age_m       | int          | YES  |     | NULL    |       |
-- | sex         | varchar(255) | YES  |     | NULL    |       |
-- | study_group | varchar(255) | YES  |     | NULL    |       |
-- +-------------+--------------+------+-----+---------+-------+
-- 7 rows in set (0.00 sec)
DESCRIBE My_samples;
-- +-----------------+--------------+------+-----+---------+-------+
-- | Field           | Type         | Null | Key | Default | Extra |
-- +-----------------+--------------+------+-----+---------+-------+
-- | manifest_id     | varchar(255) | NO   |     | NULL    |       |
-- | study_id        | varchar(255) | YES  |     | NULL    |       |
-- | visit_num       | int          | YES  |     | NULL    |       |
-- | facility        | varchar(255) | YES  |     | NULL    |       |
-- | sample_type     | varchar(255) | YES  |     | NULL    |       |
-- | tube_type       | varchar(255) | YES  |     | NULL    |       |
-- | collection_date | date         | YES  |     | NULL    |       |
-- | collection_time | time         | YES  |     | NULL    |       |
-- | aliquot_number  | int          | YES  |     | NULL    |       |
-- | aliquot_type    | varchar(255) | YES  |     | NULL    |       |
-- | aliquot_code    | varchar(255) | YES  |     | NULL    |       |
-- | initial_ml      | varchar(255) | YES  |     | NULL    |       |
-- | current_ml      | varchar(255) | YES  |     | NULL    |       |
-- | box_num         | int          | YES  |     | NULL    |       |
-- | row_num         | varchar(255) | YES  |     | NULL    |       |
-- | col_num         | int          | YES  |     | NULL    |       |
-- | hiv_status      | varchar(255) | YES  |     | NULL    |       |
-- | sample_id       | varchar(255) | NO   |     | NULL    |       |
-- | sample_name     | varchar(255) | YES  |     | NULL    |       |
-- | sample_plate    | int          | YES  |     | NULL    |       |
-- | sample_well     | int          | YES  |     | NULL    |       |
-- | i5_index_id     | varchar(255) | YES  |     | NULL    |       |
-- | index1          | varchar(255) | YES  |     | NULL    |       |
-- | i7_index_id     | varchar(255) | YES  |     | NULL    |       |
-- | index2          | varchar(255) | YES  |     | NULL    |       |
-- | sample_project  | varchar(255) | YES  |     | NULL    |       |
-- | description     | varchar(255) | YES  |     | NULL    |       |
-- +-----------------+--------------+------+-----+---------+-------+
-- 27 rows in set (0.00 sec)

-- Make a new cleaned up table from your selected columns
CREATE TABLE temp 
	SELECT My_samples.manifest_id, 
	My_samples.hiv_status, 
	My_samples.visit_num, 
	My_samples.sample_type, 
	Demographic.study_id, 
	Demographic.sex, 
	Demographic.age_y, 
	Demographic.study_group  
	FROM My_samples,Demographic 
	WHERE My_samples.study_id=Demographic.study_id;
-- Delete table
DROP TABLE temp;
-- Make a more complex table
CREATE TABLE temp 
	SELECT My_samples.manifest_id, 
	My_samples.hiv_status, 
	My_samples.visit_num, 
	My_samples.sample_type, 
	Demographic.study_id, 
	Demographic.sex, 
	Demographic.age_y, 
	Demographic.study_group,
	Diet_24hrs.fizzy_drinks,
	Diet_24hrs.antibiotics_syrup 
	FROM My_samples,Demographic,Diet_24hrs 
	WHERE My_samples.study_id=Demographic.study_id AND My_samples.study_id=Diet_24hrs.study_id;

-- write to file
-- NOTE: For now this doesn't maintain the headers, working on a way to get these included
SELECT * FROM temp ORDER BY study_id INTO OUTFILE '/Users/mann/metadata.txt' FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

-- exit mysql server
exit


