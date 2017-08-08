/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET NAMES utf8 */;
/*!50503 SET NAMES utf8mb4 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;


-- 导出 rna 的数据库结构
CREATE DATABASE IF NOT EXISTS `rna` /*!40100 DEFAULT CHARACTER SET utf8 COLLATE utf8_bin */;
USE `rna`;

-- 导出  表 rna.adar_exp 结构
CREATE TABLE IF NOT EXISTS `adar_exp` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `tissue` bigint(20) DEFAULT '0',
  `sample_type` int(11) DEFAULT '0',
  `adar1` varchar(10) COLLATE utf8_bin DEFAULT '0',
  `adar2` varchar(10) COLLATE utf8_bin DEFAULT '0',
  `editingsites` varchar(10) COLLATE utf8_bin DEFAULT '0',
  `source` varchar(10) COLLATE utf8_bin DEFAULT '0',
  KEY `rec_id` (`rec_id`),
  KEY `FK__tissue` (`tissue`),
  CONSTRAINT `FK__tissue` FOREIGN KEY (`tissue`) REFERENCES `tissue` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=159 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.candidates 结构
CREATE TABLE IF NOT EXISTS `candidates` (
  `Frequency` decimal(5,2) DEFAULT NULL,
  `Job` bigint(20) NOT NULL,
  `Position` bigint(20) NOT NULL,
  `RefSeq_feat` varchar(255) DEFAULT NULL,
  `RefSeq_gid` varchar(255) DEFAULT NULL,
  `Region` varchar(255) NOT NULL,
  `RepMask_gid` varchar(50) DEFAULT NULL,
  `Strand` varchar(1) NOT NULL,
  `ctime` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
  KEY `job` (`Job`),
  KEY `Position` (`Position`),
  KEY `Region` (`Region`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.gene_info 结构
CREATE TABLE IF NOT EXISTS `gene_info` (
  `gene_id` bigint(20) DEFAULT NULL,
  `symbol` varchar(50) DEFAULT NULL,
  `synonyms` varchar(500) DEFAULT NULL,
  `dbxrefs` text,
  `chromosome` varchar(50) DEFAULT NULL,
  `start` bigint(20) DEFAULT NULL,
  `end` bigint(20) DEFAULT NULL,
  `map_location` varchar(50) DEFAULT NULL,
  `description` text,
  `symbol_from_nomenclature_authority` varchar(50) DEFAULT NULL,
  `full_name_from_nomenclature_authority` char(255) DEFAULT NULL,
  `other_designations` text,
  KEY `gene_id` (`gene_id`),
  KEY `symbol` (`symbol`),
  KEY `synonyms` (`synonyms`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.jobs 结构
CREATE TABLE IF NOT EXISTS `jobs` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT COMMENT 'job id',
  `user` varchar(255) DEFAULT '2' COMMENT 'who submit the job',
  `st` tinyint(4) DEFAULT '0' COMMENT 'dnarna or rna denovo',
  `trace` bigint(20) DEFAULT '0' COMMENT 'id in running',
  `job` text,
  `es` text,
  `fb` tinyint(4) DEFAULT '0' COMMENT 'use hg19 rather than hg38 as ref genome',
  `status` tinyint(4) DEFAULT '0' COMMENT 'job status',
  `result` text,
  `utr` int(11) DEFAULT NULL,
  `mir` int(11) DEFAULT NULL,
  `splicing` int(11) DEFAULT NULL,
  `mis` int(11) DEFAULT NULL,
  `description` text,
  `vcode` varchar(50) DEFAULT NULL,
  `ctime` timestamp NULL DEFAULT CURRENT_TIMESTAMP,
  `adar1` varchar(10) DEFAULT '0.0',
  `adar2` varchar(10) DEFAULT '0.0',
  `tissue` varchar(10) DEFAULT NULL,
  UNIQUE KEY `id` (`id`),
  KEY `vcode` (`vcode`),
  KEY `status` (`status`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.mail_queue 结构
CREATE TABLE IF NOT EXISTS `mail_queue` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `type` tinyint(4) DEFAULT '0',
  `mail` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `vcode` varchar(50) COLLATE utf8_bin DEFAULT '0',
  `times` tinyint(4) DEFAULT '0',
  `msg` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `status` tinyint(4) DEFAULT '0',
  `ctime` datetime DEFAULT CURRENT_TIMESTAMP,
  KEY `id` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.miranda_diff 结构
CREATE TABLE IF NOT EXISTS `miranda_diff` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `sig` varchar(50) DEFAULT '0',
  `gene_symbol` varchar(50) DEFAULT '0',
  KEY `sig` (`sig`),
  KEY `gene_symbol` (`gene_symbol`),
  KEY `rec_id` (`rec_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.miranda_score 结构
CREATE TABLE IF NOT EXISTS `miranda_score` (
  `rec_id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `mirna` varchar(255) COLLATE utf8_bin NOT NULL,
  `sig` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `job` bigint(20) NOT NULL,
  `tag` varchar(255) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `way` varchar(10) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `gene_symbol` varchar(255) COLLATE utf8_bin NOT NULL,
  `score` decimal(5,2) NOT NULL,
  `energy` decimal(5,2) NOT NULL,
  `mi_start` int(11) NOT NULL,
  `mi_end` int(11) NOT NULL,
  `utr_start` int(11) NOT NULL,
  `utr_end` int(11) NOT NULL,
  `match_len` int(11) NOT NULL,
  `identity` varchar(10) COLLATE utf8_bin NOT NULL,
  `similarity` varchar(10) COLLATE utf8_bin NOT NULL,
  `mir_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  `lines` varchar(255) COLLATE utf8_bin NOT NULL,
  `utr_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  UNIQUE KEY `rec_id` (`rec_id`),
  KEY `mirna` (`mirna`),
  KEY `tag` (`tag`),
  KEY `gene_symbol` (`gene_symbol`),
  KEY `sig` (`sig`),
  KEY `job` (`job`),
  KEY `wag` (`way`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.mirediting 结构
CREATE TABLE IF NOT EXISTS `mirediting` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `mirna` varchar(255) DEFAULT '0',
  `accession` varchar(255) DEFAULT '0',
  `edit_pos_raw` varchar(255) DEFAULT '0',
  `edit_pos_chr` varchar(255) DEFAULT '0',
  `event` varchar(255) DEFAULT '0',
  `role` varchar(255) DEFAULT '0',
  `sequence` varchar(255) DEFAULT '0',
  `old_t` int(11) DEFAULT '0',
  `new_t` int(11) DEFAULT '0',
  `com_t` int(11) DEFAULT '0',
  `sig` varchar(255) DEFAULT '0',
  `job` bigint(20) DEFAULT '0',
  KEY `sig` (`sig`),
  KEY `job` (`job`),
  KEY `Column 1` (`rec_id`)
) ENGINE=InnoDB AUTO_INCREMENT=241 DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.mirgenome 结构
CREATE TABLE IF NOT EXISTS `mirgenome` (
  `chromosome` varchar(2) COLLATE utf8_bin NOT NULL,
  `source` varchar(50) COLLATE utf8_bin NOT NULL,
  `type` varchar(50) COLLATE utf8_bin NOT NULL,
  `start` bigint(20) NOT NULL,
  `end` bigint(20) NOT NULL,
  `score` varchar(50) COLLATE utf8_bin NOT NULL,
  `strand` varchar(1) COLLATE utf8_bin NOT NULL,
  `phase` varchar(50) COLLATE utf8_bin NOT NULL,
  `mirna_id` varchar(50) COLLATE utf8_bin NOT NULL,
  `alias` varchar(50) COLLATE utf8_bin NOT NULL,
  `mirna_name` varchar(50) COLLATE utf8_bin NOT NULL,
  `derives_from` varchar(50) COLLATE utf8_bin DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.mirmap_score 结构
CREATE TABLE IF NOT EXISTS `mirmap_score` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `mirna` varchar(50) DEFAULT '0',
  `sig` varchar(50) DEFAULT '0',
  `job` varchar(50) DEFAULT '0',
  `tag` varchar(50) DEFAULT '0',
  `way` varchar(50) DEFAULT '0',
  `gene_symbol` varchar(50) DEFAULT '0',
  `transcript_id` varchar(50) DEFAULT '0',
  `dg_duplex` varchar(50) DEFAULT '0',
  `dg_binding` varchar(50) DEFAULT '0',
  `dg_duplex_seed` varchar(50) DEFAULT '0',
  `dg_binding_seed` varchar(50) DEFAULT '0',
  `utr_start` varchar(50) DEFAULT '0',
  `utr_end` varchar(50) DEFAULT '0',
  `utr3` varchar(50) DEFAULT '0',
  KEY `rec_id` (`rec_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.mirna 结构
CREATE TABLE IF NOT EXISTS `mirna` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `accession` varchar(50) COLLATE utf8_bin NOT NULL,
  `mir_id` varchar(50) COLLATE utf8_bin NOT NULL,
  `status` varchar(50) COLLATE utf8_bin NOT NULL,
  `sequence` text COLLATE utf8_bin NOT NULL,
  `mature1_acc` varchar(50) COLLATE utf8_bin NOT NULL,
  `mature1_id` varchar(50) COLLATE utf8_bin NOT NULL,
  `mature1_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  `mature2_acc` varchar(50) COLLATE utf8_bin NOT NULL,
  `mature2_id` varchar(50) COLLATE utf8_bin NOT NULL,
  `mature2_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  `symbol` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  KEY `rec_id` (`rec_id`),
  KEY `mature1_acc` (`mature1_acc`),
  KEY `mature1_id` (`mature1_id`),
  KEY `mature2_acc` (`mature2_acc`),
  KEY `mature2_id` (`mature2_id`),
  KEY `mir_id` (`mir_id`)
) ENGINE=InnoDB AUTO_INCREMENT=28343 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.mirnatarget 结构
CREATE TABLE IF NOT EXISTS `mirnatarget` (
  `mirna_id` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `mature_name` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `transcript_id` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `transcript_stable_id` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `transcript_canonical` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `transcript_chr` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `transcript_strand` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `gene_stable_id` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `gene_name` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `name` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `site_id` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `site_end` bigint(20) DEFAULT NULL,
  `site_chr_end` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `seed_length` int(11) DEFAULT NULL,
  `seed_start` bigint(20) DEFAULT NULL,
  `seed_end` bigint(20) DEFAULT NULL,
  `seed_mismatches_nogu` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `seed_gu` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `tgs_au` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `tgs_position` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `tgs_pairing3p` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `dg_duplex` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `dg_binding` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `dg_duplex_seed` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `dg_bingding_seed` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `dg_open` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `dg_total` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `prob_exact` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `prob_binomial` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `cons_bls` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `slec_phylop` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `mirmap_score` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  KEY `mature_name` (`mature_name`),
  KEY `gene_name` (`gene_name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.mirna_map 结构
CREATE TABLE IF NOT EXISTS `mirna_map` (
  `acc` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `mir` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  KEY `acc` (`acc`),
  KEY `mir` (`mir`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.mir_overview 结构
CREATE TABLE IF NOT EXISTS `mir_overview` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `job_id` bigint(20) NOT NULL DEFAULT '0',
  `reg` longtext NOT NULL,
  `rel` longtext NOT NULL,
  KEY `rec_id` (`rec_id`),
  KEY `job_id` (`job_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.missense 结构
CREATE TABLE IF NOT EXISTS `missense` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `job` bigint(20) NOT NULL DEFAULT '0',
  `chromosome` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `position` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `gene` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `transcript` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `relpos` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `fr` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `to` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `ctime` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`rec_id`),
  KEY `job` (`job`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.searchable 结构
CREATE TABLE IF NOT EXISTS `searchable` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `sig` varchar(50) DEFAULT NULL,
  `Frequency` decimal(5,2) DEFAULT NULL,
  `Job` bigint(20) NOT NULL,
  `Position` bigint(20) NOT NULL,
  `RefSeq_feat` varchar(255) DEFAULT NULL,
  `RefSeq_gid` varchar(255) DEFAULT NULL,
  `Region` varchar(255) NOT NULL,
  `RepMask_gid` varchar(50) DEFAULT NULL,
  `Strand` varchar(1) NOT NULL,
  `ctime` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `job` (`Job`),
  KEY `Position` (`Position`),
  KEY `Region` (`Region`),
  KEY `RefSeq_gid` (`RefSeq_gid`),
  KEY `列 2` (`sig`)
) ENGINE=InnoDB AUTO_INCREMENT=1730404 DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.smir 结构
CREATE TABLE IF NOT EXISTS `smir` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `mirna` varchar(255) DEFAULT '0',
  `accession` varchar(255) DEFAULT '0',
  `edit_pos_raw` varchar(255) DEFAULT '0',
  `edit_pos_chr` varchar(255) DEFAULT '0',
  `event` varchar(255) DEFAULT '0',
  `role` varchar(255) DEFAULT '0',
  `sequence` varchar(255) DEFAULT '0',
  `old_t` int(11) DEFAULT '0',
  `new_t` int(11) DEFAULT '0',
  `com_t` int(11) DEFAULT '0',
  `sig` varchar(255) DEFAULT '0',
  `job` bigint(20) DEFAULT '1',
  UNIQUE KEY `Column 1` (`rec_id`),
  KEY `job` (`job`),
  KEY `edit_pos_raw` (`edit_pos_raw`)
) ENGINE=InnoDB AUTO_INCREMENT=1039 DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.smiranda 结构
CREATE TABLE IF NOT EXISTS `smiranda` (
  `rec_id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `mirna` varchar(255) COLLATE utf8_bin NOT NULL,
  `sig` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `job` bigint(20) NOT NULL DEFAULT '1',
  `tag` varchar(255) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `way` varchar(10) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `gene_symbol` varchar(255) COLLATE utf8_bin NOT NULL,
  `score` decimal(5,2) NOT NULL,
  `energy` decimal(5,2) NOT NULL,
  `mi_start` int(11) NOT NULL,
  `mi_end` int(11) NOT NULL,
  `utr_start` int(11) NOT NULL,
  `utr_end` int(11) NOT NULL,
  `match_len` int(11) NOT NULL,
  `identity` varchar(10) COLLATE utf8_bin NOT NULL,
  `similarity` varchar(10) COLLATE utf8_bin NOT NULL,
  `mir_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  `lines` varchar(255) COLLATE utf8_bin NOT NULL,
  `utr_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  UNIQUE KEY `rec_id` (`rec_id`),
  KEY `mirna` (`mirna`),
  KEY `tag` (`tag`),
  KEY `gene_symbol` (`gene_symbol`),
  KEY `sig` (`sig`),
  KEY `job` (`job`),
  KEY `wag` (`way`)
) ENGINE=InnoDB AUTO_INCREMENT=426293 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.smis 结构
CREATE TABLE IF NOT EXISTS `smis` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `job` bigint(20) NOT NULL DEFAULT '0',
  `chromosome` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `position` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `gene` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `transcript` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `relpos` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `fr` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `to` varchar(50) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `ctime` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`rec_id`),
  KEY `job` (`job`)
) ENGINE=InnoDB AUTO_INCREMENT=3875 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.sms 结构
CREATE TABLE IF NOT EXISTS `sms` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `mirna` varchar(50) DEFAULT '0',
  `sig` varchar(50) DEFAULT '0',
  `job` varchar(50) DEFAULT '0',
  `tag` varchar(50) DEFAULT '0',
  `way` varchar(50) DEFAULT '0',
  `gene_symbol` varchar(50) DEFAULT '0',
  `transcript_id` varchar(50) DEFAULT '0',
  `dg_duplex` varchar(50) DEFAULT '0',
  `dg_binding` varchar(50) DEFAULT '0',
  `dg_duplex_seed` varchar(50) DEFAULT '0',
  `dg_binding_seed` varchar(50) DEFAULT '0',
  `utr_start` varchar(50) DEFAULT '0',
  `utr_end` varchar(50) DEFAULT '0',
  `utr3` varchar(50) DEFAULT '0',
  KEY `rec_id` (`rec_id`)
) ENGINE=InnoDB AUTO_INCREMENT=1255632 DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.splicing_event 结构
CREATE TABLE IF NOT EXISTS `splicing_event` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `gene` varchar(255) DEFAULT '0',
  `chromosome` varchar(255) DEFAULT '0',
  `pos` bigint(20) DEFAULT '0',
  `type` tinyint(4) DEFAULT '0',
  `raw_score` decimal(6,3) DEFAULT '0.000',
  `new_score` decimal(6,3) DEFAULT '0.000',
  `variation` decimal(6,3) DEFAULT '0.000',
  `order` tinyint(4) DEFAULT '0',
  `job` bigint(20) DEFAULT '0',
  `transcript` varchar(50) DEFAULT NULL,
  KEY `rec_id` (`rec_id`),
  KEY `gene` (`gene`),
  KEY `pos` (`pos`),
  KEY `chromosome` (`chromosome`)
) ENGINE=InnoDB AUTO_INCREMENT=7471 DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.sse 结构
CREATE TABLE IF NOT EXISTS `sse` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `gene` varchar(255) DEFAULT '0',
  `pos` bigint(20) DEFAULT '0',
  `type` tinyint(4) DEFAULT '0',
  `raw_score` decimal(6,3) DEFAULT '0.000',
  `new_score` decimal(6,3) DEFAULT '0.000',
  `variation` decimal(6,3) DEFAULT '0.000',
  `order` tinyint(4) DEFAULT '0',
  `job` bigint(20) DEFAULT '0',
  `transcript` varchar(50) DEFAULT NULL,
  KEY `rec_id` (`rec_id`),
  KEY `job` (`job`),
  KEY `pos` (`pos`),
  KEY `gene` (`gene`)
) ENGINE=InnoDB AUTO_INCREMENT=269653 DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.sum 结构
CREATE TABLE IF NOT EXISTS `sum` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `gene` varchar(50) DEFAULT '0',
  `sig` varchar(50) DEFAULT '0',
  `job` varchar(50) NOT NULL DEFAULT '0',
  `tag` varchar(50) DEFAULT '0',
  `way` varchar(50) DEFAULT '0',
  `mir` varchar(50) DEFAULT '0',
  `dg_duplex` varchar(50) DEFAULT '0',
  `dg_binding` varchar(50) DEFAULT '0',
  `dg_duplex_seed` varchar(50) DEFAULT '0',
  `dg_binding_seed` varchar(50) DEFAULT '0',
  `utr_start` varchar(50) DEFAULT '0',
  `utr_end` varchar(50) DEFAULT '0',
  `utr3` varchar(50) DEFAULT '0',
  KEY `rec_id` (`rec_id`),
  KEY `sig` (`sig`),
  KEY `tag` (`tag`),
  KEY `job` (`job`)
) ENGINE=InnoDB AUTO_INCREMENT=224363623 DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
-- 导出  表 rna.sutr3 结构
CREATE TABLE IF NOT EXISTS `sutr3` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `gene` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `edit_pos_chr` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `edit_pos_raw` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `event` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `old_t` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `new_t` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `com_t` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `sig` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `job` varchar(255) COLLATE utf8_bin DEFAULT '0',
  UNIQUE KEY `rec_id` (`rec_id`),
  KEY `job` (`job`),
  KEY `edit_pos_raw` (`edit_pos_raw`),
  KEY `gene` (`gene`),
  KEY `edit_pos_chr` (`edit_pos_chr`)
) ENGINE=InnoDB AUTO_INCREMENT=110202 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.tissue 结构
CREATE TABLE IF NOT EXISTS `tissue` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `tissue` varchar(50) COLLATE utf8_bin DEFAULT '0',
  `status` tinyint(4) DEFAULT '1',
  UNIQUE KEY `id` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=32 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.utr3 结构
CREATE TABLE IF NOT EXISTS `utr3` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `chr` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `start` bigint(20) DEFAULT NULL,
  `end` bigint(20) DEFAULT NULL,
  `infasta` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `strand` varchar(1) COLLATE utf8_bin DEFAULT NULL,
  `gene_symbol` varchar(255) COLLATE utf8_bin DEFAULT NULL,
  KEY `rec_id` (`rec_id`),
  KEY `chr` (`chr`),
  KEY `start` (`start`),
  KEY `end` (`end`),
  KEY `strand` (`strand`),
  KEY `infasta` (`infasta`)
) ENGINE=InnoDB AUTO_INCREMENT=77078 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.utr3editing 结构
CREATE TABLE IF NOT EXISTS `utr3editing` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `gene` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `edit_pos_chr` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `edit_pos_raw` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `event` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `old_t` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `new_t` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `com_t` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `sig` varchar(255) COLLATE utf8_bin DEFAULT '0',
  `job` varchar(255) COLLATE utf8_bin DEFAULT '0',
  UNIQUE KEY `rec_id` (`rec_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.utr3_prediction 结构
CREATE TABLE IF NOT EXISTS `utr3_prediction` (
  `rec_id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `mirna` varchar(255) COLLATE utf8_bin NOT NULL,
  `fastatag` varchar(255) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `repeat_chk` varchar(255) COLLATE utf8_bin NOT NULL,
  `sig` varchar(50) COLLATE utf8_bin DEFAULT NULL,
  `tag` varchar(255) COLLATE utf8_bin NOT NULL DEFAULT '0',
  `gene_symbol` varchar(255) COLLATE utf8_bin NOT NULL,
  `score` decimal(5,2) NOT NULL,
  `energy` decimal(5,2) NOT NULL,
  `mi_start` int(11) NOT NULL,
  `mi_end` int(11) NOT NULL,
  `utr_start` int(11) NOT NULL,
  `utr_end` int(11) NOT NULL,
  `match_len` int(11) NOT NULL,
  `identity` varchar(50) COLLATE utf8_bin NOT NULL,
  `similarity` varchar(50) COLLATE utf8_bin NOT NULL,
  `mir_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  `lines` varchar(255) COLLATE utf8_bin NOT NULL,
  `utr_seq` varchar(255) COLLATE utf8_bin NOT NULL,
  `trace` bigint(20) unsigned NOT NULL DEFAULT '0',
  UNIQUE KEY `rec_id` (`rec_id`),
  KEY `mirna` (`mirna`),
  KEY `tag` (`tag`),
  KEY `gene_symbol` (`gene_symbol`),
  KEY `sig` (`sig`),
  KEY `repeat_chk` (`repeat_chk`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- 数据导出被取消选择。
-- 导出  表 rna.utr_map 结构
CREATE TABLE IF NOT EXISTS `utr_map` (
  `rec_id` bigint(20) NOT NULL AUTO_INCREMENT,
  `gene` varchar(50) DEFAULT '0',
  `sig` varchar(50) DEFAULT '0',
  `job` varchar(50) NOT NULL DEFAULT '0',
  `tag` varchar(50) DEFAULT '0',
  `way` varchar(50) DEFAULT '0',
  `mir` varchar(50) DEFAULT '0',
  `dg_duplex` varchar(50) DEFAULT '0',
  `dg_binding` varchar(50) DEFAULT '0',
  `dg_duplex_seed` varchar(50) DEFAULT '0',
  `dg_binding_seed` varchar(50) DEFAULT '0',
  `utr_start` varchar(50) DEFAULT '0',
  `utr_end` varchar(50) DEFAULT '0',
  `utr3` varchar(50) DEFAULT '0',
  KEY `rec_id` (`rec_id`),
  KEY `sig` (`sig`),
  KEY `tag` (`tag`),
  KEY `job` (`job`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- 数据导出被取消选择。
/*!40101 SET SQL_MODE=IFNULL(@OLD_SQL_MODE, '') */;
/*!40014 SET FOREIGN_KEY_CHECKS=IF(@OLD_FOREIGN_KEY_CHECKS IS NULL, 1, @OLD_FOREIGN_KEY_CHECKS) */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
