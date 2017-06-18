<?php
namespace Viewer\Controller;
use Think\Controller;
class Utr3Controller extends Controller {
    public function index(){
        
    }
	
	public function pres($sig, $sf='n'){
		if ($sf == 'y'){
			$overDB = M('sutr3');
			$db_targets = M('sum');
			$dbs = M('searchable');
			session("search", 1);
		}else{
			$overDB = M('utr3editing');
			$db_targets = M('utr_map');
			$dbs = M('candidates');
			session("search", 0);
		}
		if (session('jobId')) {
			//$overDB = M('utr3editing');
			//$overDB = M('sutr3');
			$map = array('sig'	=>	$sig);
			$over = $overDB->where($map)->find();
			if ($over) {
				$bi = $dbs->where("job=".session('jobId')." AND position=".$over['edit_pos_raw']." AND region='".$over['edit_pos_chr']."'")->find();
				#$Model = new \Think\Model();
				//$db_targets = M('utr_map');
				#$db_targets = M('sum');
				$oldTargets = $db_targets->field('rec_id, mir')->where(array('sig'=>$sig,'tag'=>'RAW','job'=>session('jobId')))->order('mir desc')->limit(100)->select();
				$newTargets = $db_targets->field('rec_id, mir')->where(array('sig'=>$sig,'tag'=>'NEW','job'=>session('jobId')))->order('mir desc')->limit(100)->select();
				#$newTargets = $Model->query("SELECT rec_id, mir FROM utr_map WHERE tag = 'NEW' AND sig = '{$sig}' ORDER BY score desc LIMIT 100;");
				$inactiveTargets = $db_targets->field('rec_id, mir')->where(array('sig'=>$sig,'tag'=>'DIE','job'=>session('jobId')))->order('mir desc')->limit(100)->select();
				#$inactiveTargets = $Model->query("SELECT rec_id, mirna FROM utr3_prediction WHERE tag = 'RAW' AND sig = '{$sig}' AND mirna NOT IN (SELECT mirna FROM utr3_prediction WHERE tag != 'RAW' AND sig = '{$sig}' GROUP BY mirna) GROUP BY mirna ORDER BY score desc LIMIT 100;");
				$this->assign('hit', $over);
				$this->assign('bi', $bi);
				$this->assign("oldList", $oldTargets);
				$this->assign("newList", $newTargets);
				$this->assign("incList", $inactiveTargets);
				$this->display();
			}else{
				$this->error("No sites found!");
			}
		}else{
			$this->error("> Are you attacking RNA Editing Plus server?<br />This form is under Bro Li's protection.");
		}
	}
	
	public function pair($rec){
		if (session("search") == 1){
			$db = M('sum');
		}else{
			$db = M('utr_map');
		}
		
		$ret = $db->where(array('rec_id'=>intval($rec)))->select();
		if ($ret){
			$this->assign("hit", $ret[0]);
			$this->display();
		}else{
			$this->error("No record found!");
		}
	}
	
	public function priOverview($rec){
		if (session('jobId')) {
			$geneInfo = M('gene_info');
			$map = array('symbol'	=>	$rec);
			$gi = $geneInfo->where($map)->find();
			if (session("jobId") == 1){
				$miEd = M('sutr3');
			}else{
				$miEd = M('utr3editing');
			}
			$emap = array(
					'gene'	=>	$rec, 
					'job'	=>	session('jobId')
			);
			$ei = $miEd->where($emap)->select();
			$this->assign('gi', $gi);
			$this->assign('rlist', $ei);
			$this->display();
		}else{
			$this->error("> Are you attacking RNA Editing Plus server?<br />This form is under Bro Li's protection.");
		}
	}
	
	public function geneList($sig, $type, $p=1){
		if (session("search") == 1){
			$db = M('sum');
		}else{
			$db = M('utr_map');
		}

		//$db = M('utr_map');
		#$db = M('sum');
		$jobId = session('jobId');
		if ($jobId){
			if ($type == 'raw') {
				$map = array('sig'=>$sig, 'tag'=>'RAW');
			}elseif($type == 'mut'){
				$map = array('sig'=>$sig, 'tag'=>'NEW');
			}elseif($type == 'inc'){
				$map = array('sig'=>$sig, 'tag'=>'DIE');
			}else{
				$this->error('> Are you attacking RNA Editing Plus server?');
			}
			$ret = $db->where($map)->order('mir desc')->page($_GET['p'].',20')->select();
			$count      = $db->where($map)->count();// 查询满足要求的总记录数
			$Page       = new \Think\Page($count,20);// 实例化分页类 传入总记录数和每页显示的记录数
			$show       = $Page->show();// 分页显示输出
			$this->assign('page', $show);
			$this->assign("rlist", $ret);
			$this->assign("sig", $sig);
			$this->assign("gene", $ret[0]["gene"]);
			$this->display("pairList");
		}else{
			$this->error("> Are you attacking RNA Editing Plus server?<br />This form is under Bro Li's protection.");
		}
		
	}
	
	public function search($symbol, $type, $sig){
		$symbol = str_replace('hsa-', '', $symbol);
		if (session("search") == 1){
			$db = M('sum');
		}else{
			$db = M('utr_map');
		}
		switch ($type)
		{
			case 'RAW':
				$map = array(
						'mir'	=>	$symbol,
						'sig'			=>	$sig,
						'tag'			=>	'RAW',
				);
				$ret = $db->where($map)->find();
				break;
			case 'MUT':
				/*
				$Model = new \Think\Model();
				$ret = $Model->query("SELECT * FROM utr3_prediction WHERE tag != 'RAW' AND sig = '{$sig}' AND mirna = '{$symbol}' AND mirna NOT IN (SELECT mirna FROM utr3_prediction WHERE tag = 'RAW' AND sig = '{$sig}' GROUP BY mirna);");
				*/
				$map = array(
						'mir'	=>	$symbol,
						'sig'	=>	$sig,
						'tag'	=>	'NEW',
				);
				$ret = $db->where($map)->find();
				break;
			case 'INC':
				/*
				$Model = new \Think\Model();
				$ret = $Model->query("SELECT * FROM utr3_prediction WHERE tag = 'RAW' AND sig = '{$sig}' AND mirna = '{$symbol}' AND mirna NOT IN (SELECT mirna FROM utr3_prediction WHERE tag != 'RAW' AND sig = '{$sig}' GROUP BY mirna);");
				*/
				$map = array(
						'mir'	=>	$symbol,
						'sig'	=>	$sig,
						'tag'	=>	'DIE',
				);
				$ret = $db->where($map)->find();
				break;
			default:
		}
		if ($ret){
			$this->assign("hit", $ret);
			echo $this->fetch();
			//var_dump ($this->fetch());
		}else{
			echo "No record found!";
		}
	
	}
}
?>