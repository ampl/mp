-- MySQL dump 10.13  Distrib 5.5.28, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: test
-- ------------------------------------------------------
-- Server version	5.5.28-0ubuntu0.12.10.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `Amounts`
--

DROP TABLE IF EXISTS `Amounts`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Amounts` (
  `NUTR` varchar(4) DEFAULT NULL,
  `FOOD` varchar(5) DEFAULT NULL,
  `amt` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Amounts`
--

LOCK TABLES `Amounts` WRITE;
/*!40000 ALTER TABLE `Amounts` DISABLE KEYS */;
INSERT INTO `Amounts` VALUES ('A','BEEF',60),('C','BEEF',20),('B1','BEEF',10),('B2','BEEF',15),('NA','BEEF',938),('CAL','BEEF',295),('A','CHK',8),('C','CHK',0),('B1','CHK',20),('B2','CHK',20),('NA','CHK',2180),('CAL','CHK',770),('A','FISH',8),('C','FISH',10),('B1','FISH',15),('B2','FISH',10),('NA','FISH',945),('CAL','FISH',440),('A','HAM',40),('C','HAM',40),('B1','HAM',35),('B2','HAM',10),('NA','HAM',278),('CAL','HAM',430),('A','MCH',15),('C','MCH',35),('B1','MCH',15),('B2','MCH',15),('NA','MCH',1182),('CAL','MCH',315),('A','MTL',70),('C','MTL',30),('B1','MTL',15),('B2','MTL',15),('NA','MTL',896),('CAL','MTL',400),('A','SPG',25),('C','SPG',50),('B1','SPG',25),('B2','SPG',15),('NA','SPG',1329),('CAL','SPG',370),('A','TUR',60),('C','TUR',20),('B1','TUR',15),('B2','TUR',10),('NA','TUR',1397),('CAL','TUR',450);
/*!40000 ALTER TABLE `Amounts` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Foods`
--

DROP TABLE IF EXISTS `Foods`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Foods` (
  `FOOD` varchar(5) DEFAULT NULL,
  `cost` double DEFAULT NULL,
  `f_min` double DEFAULT NULL,
  `f_max` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Foods`
--

LOCK TABLES `Foods` WRITE;
/*!40000 ALTER TABLE `Foods` DISABLE KEYS */;
INSERT INTO `Foods` VALUES ('BEEF',3.19,2,10),('CHK',2.59,2,10),('FISH',2.29,2,10),('HAM',2.89,2,10),('MCH',1.89,2,10),('MTL',1.99,2,10),('SPG',1.99,2,10),('TUR',2.49,2,10);
/*!40000 ALTER TABLE `Foods` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Nutrients`
--

DROP TABLE IF EXISTS `Nutrients`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Nutrients` (
  `NUTR` varchar(4) DEFAULT NULL,
  `n_min` double DEFAULT NULL,
  `n_max` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Nutrients`
--

LOCK TABLES `Nutrients` WRITE;
/*!40000 ALTER TABLE `Nutrients` DISABLE KEYS */;
INSERT INTO `Nutrients` VALUES ('A',700,20000),('C',700,20000),('B1',700,20000),('B2',700,20000),('NA',0,50000),('CAL',16000,24000);
/*!40000 ALTER TABLE `Nutrients` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2012-11-06  9:03:46
