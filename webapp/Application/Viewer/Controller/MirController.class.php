<?php
namespace Viewer\Controller;
use Think\Controller;
class MirController extends Controller {
    public function index(){
        
    }
	
	public function pres($rec, $sf='n'){
		if ($sf == 'y'){
			$db = M('smir');
			$db_targets_seed = M('sms');
			$dbs = M('searchable');
			$db_targets_mature = M('smiranda');
			session("search", 1);
		}else{
			$db = M('mirediting');
			$db_targets_seed = M('mirmap_score');
			$dbs = M('candidates');
			$db_targets_mature = M('miranda_score');
			session("search", 0);
		}
		
		$ret = $db->where(array('rec_id'=>intval($rec), 'job'	=>	session('jobId')))->select();

		$sig = $ret[0]['sig'];
		if ($ret[0]['role'] != 0){
			$oldTargets = $db_targets_seed->field('rec_id, gene_symbol')->where(array('sig'=>$sig,'tag'=>'RAW'))->order('gene_symbol asc')->limit(100)->group('gene_symbol')->select();
			//$Model = new \Think\Model();
			$newTargets = $db_targets_seed->field('rec_id, gene_symbol')->where(array('sig'=>$sig,'tag'=>'NEW'))->order('gene_symbol asc')->limit(100)->group('gene_symbol')->select();
			#$newTargets = $Model->query("SELECT rec_id, gene_symbol FROM mirmap_score WHERE tag != 'RAW' AND AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM mirmap_score WHERE tag = 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY gene_symbol asc LIMIT 100;");
			#$inactiveTargets = $Model->query("SELECT rec_id, gene_symbol FROM mirmap_score WHERE tag = 'RAW' AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM mirmap_score WHERE tag != 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY gene_symbol asc LIMIT 100;");
			$inactiveTargets = $db_targets_seed->field('rec_id, gene_symbol')->where(array('sig'=>$sig,'tag'=>'DIE'))->order('gene_symbol asc')->limit(100)->group('gene_symbol')->select();
			if ($ret){
				$bi = $dbs->where("job=".session('jobId')." AND position=".$ret[0]['edit_pos_raw']." AND region='".$ret[0]['edit_pos_chr']."'")->find();
				$this->assign("hit", $ret[0]);
				$this->assign('bi', $bi);
				$this->assign("oldList", $oldTargets);
				$this->assign("newList", $newTargets);
				$this->assign("incList", $inactiveTargets);
				$this->display('pres_seed');
			}else{
				$this->error("No sites found!");
			}
		}else{
			$oldTargets = $db_targets_mature->field('rec_id, gene_symbol')->where(array('sig'=>$sig,'tag'=>'RAW'))->order('score desc')->limit(100)->group('gene_symbol')->select();
			$newTargets = $db_targets_mature->field('rec_id, gene_symbol')->where(array('sig'=>$sig,'way'=>'NEW'))->order('score desc')->limit(100)->group('gene_symbol')->select();
			$inactiveTargets = $db_targets_mature->field('rec_id, gene_symbol')->where(array('sig'=>$sig,'way'=>'DIE'))->order('score desc')->limit(100)->group('gene_symbol')->select();
			//$Model = new \Think\Model();
			//$newTargets = $Model->query("SELECT rec_id, gene_symbol FROM miranda_score WHERE tag != 'RAW' AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM miranda_score WHERE tag = 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY score desc LIMIT 100;");
			//$inactiveTargets = $Model->query("SELECT rec_id, gene_symbol FROM miranda_score WHERE tag = 'RAW' AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM miranda_score WHERE tag != 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY score desc LIMIT 100;");
			if ($ret){
				$bi = $dbs->where("job=".session('jobId')." AND position=".$ret[0]['edit_pos_raw']." AND region='".$ret[0]['edit_pos_chr']."'")->find();
				$this->assign("hit", $ret[0]);
				$this->assign('bi', $bi);
				$this->assign("oldList", $oldTargets);
				$this->assign("newList", $newTargets);
				$this->assign("incList", $inactiveTargets);
				$this->display();
			}else{
				$this->error("No sites found!");
			}
		}
	}
	public function priOverview($rec){
		if (session('jobId')) {
			$mirna = M('mirgenome');
			$map = array('mirna_name'	=>	$rec);
			$gi = $mirna->where($map)->find();
			
			if (session("search") == 1){
				$miEd = M('smir');
			}else{
				$miEd = M('mirediting');
			}
			$emap = array('mirna'	=>	$rec, 'job'	=>	session('jobId'));
			$ei = $miEd->where($emap)->select();
			$this->assign('gi', $gi);
			$this->assign('rlist', $ei);
			$this->display();
		}else{
			$this->error("> Are you attacking RNA Editing Plus server?<br />This form is under Bro Li's protection.");
		}
	}
	public function priPres($rec){
		if (session("search") == 1){
			$db = M('smir');
			$db_targets = M('smiranda');
		}else{
			$db = M('mirediting');
			$db_targets = M('miranda_score');
		}
		$ret = $db->where(array('mirna'=>$rec, 'job'=>session('jobId')))->select();
		$sig = $ret[0]['sig'];
		$oldTargets = $db_targets->field('rec_id, gene_symbol')->where(array('sig'=>$sig,'tag'=>'RAW'))->order('score desc')->limit(100)->group('gene_symbol')->select();
		$Model = new \Think\Model();
		$newTargets = $Model->query("SELECT rec_id, gene_symbol FROM miranda_score WHERE tag != 'RAW' AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM miranda_score WHERE tag = 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY score desc LIMIT 100;");
		$inactiveTargets = $Model->query("SELECT rec_id, gene_symbol FROM miranda_score WHERE tag = 'RAW' AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM miranda_score WHERE tag != 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY score desc LIMIT 100;");
		if ($ret){
			$this->assign("hit", $ret[0]);
			$this->assign("oldList", $oldTargets);
			$this->assign("newList", $newTargets);
			$this->assign("incList", $inactiveTargets);
			$this->display('pres');
		}else{
			$this->error("No sites found!");
		}
	}
	
	public function pair($rec){
		if (session("search") == 1){
			$db = M('smiranda');
		}else{
			$db = M('miranda_score');
		}
		
		$ret = $db->where(array('rec_id'=>intval($rec)))->select();
		if ($ret){
			$this->assign("hit", $ret[0]);
			$this->display();
		}else{
			$this->error("No record found!");
		}
	}
	
	public function pairS($rec){
		if (session("search") == 1){
			$db = M('sms');
		}else{
			$db = M('mirmap_score');
		}
		$ret = $db->where(array('rec_id'=>intval($rec)))->select();
		if ($ret){
			$this->assign("hit", $ret[0]);
			$this->display();
		}else{
			$this->error("No record found!");
		}
	}
	
	public function search($symbol, $type, $sig){
		if (session("search") == 1){
			$db = M('smiranda');
		}else{
			$db = M('miranda_score');
		}
		switch ($type)
		{
			case 'RAW':
				#$db = M('mirmap_score');
				$map = array(
							'gene_symbol'	=>	$symbol,
							'sig'			=>	$sig,
							'tag'			=>	'RAW',
				);
				$ret = $db->where($map)->select();
				break;
			case 'MUT':
				$map = array(
							'gene_symbol'	=>	$symbol,
							'sig'			=>	$sig,
							'tag'			=>	'NEW',
				);
				$ret = $db->where($map)->select();
				#$Model = new \Think\Model();
				#$ret = $Model->query("SELECT * FROM mirmap_score WHERE tag = 'NEW' AND sig = '{$sig}' AND gene_symbol = '{$symbol}';");
				break;
			case 'INC':
				$map = array(
							'gene_symbol'	=>	$symbol,
							'sig'			=>	$sig,
							'tag'			=>	'DIE',
				);
				$ret = $db->where($map)->select();
				#$Model = new \Think\Model();
				#$ret = $Model->query("SELECT * FROM mirmap_score WHERE tag = 'DIE' AND sig = '{$sig}' AND gene_symbol = '{$symbol}';");
				break;
			default:
		}
		if ($ret){
			$this->assign("hitList", $ret);
			echo $this->fetch();
		}else{
			echo "No record found!";
		}
	}
	public function searchS($symbol, $type, $sig){
		if (session("search") == 1){
			$db = M('sms');
		}else{
			$db = M('mirmap_score');
		}
		switch ($type)
		{
			case 'RAW':
				#$db = M('mirmap_score');
				$map = array(
							'gene_symbol'	=>	$symbol,
							'sig'			=>	$sig,
							'tag'			=>	'RAW',
				);
				$ret = $db->where($map)->select();
				break;
			case 'MUT':
				$map = array(
							'gene_symbol'	=>	$symbol,
							'sig'			=>	$sig,
							'tag'			=>	'NEW',
				);
				$ret = $db->where($map)->select();
				#$Model = new \Think\Model();
				#$ret = $Model->query("SELECT * FROM mirmap_score WHERE tag = 'NEW' AND sig = '{$sig}' AND gene_symbol = '{$symbol}';");
				break;
			case 'INC':
				$map = array(
							'gene_symbol'	=>	$symbol,
							'sig'			=>	$sig,
							'tag'			=>	'DIE',
				);
				$ret = $db->where($map)->select();
				#$Model = new \Think\Model();
				#$ret = $Model->query("SELECT * FROM mirmap_score WHERE tag = 'DIE' AND sig = '{$sig}' AND gene_symbol = '{$symbol}';");
				break;
			default:
		}
		if ($ret){
			$this->assign("hitList", $ret);
			echo $this->fetch();
		}else{
			echo "No record found!";
		}
	}
	public function geneList($sig, $type){
		if (session("search") == 1){
			$db = M('smiranda');
		}else{
			$db = M('miranda_score');
		}
		if ($type == 'raw') {
			$ret = $db->where(array('sig'=>$sig, 'tag'=>'RAW'))->select();
		}elseif($type == 'mut'){
			$Model = new \Think\Model();
			$ret = $Model->query("SELECT * FROM miranda_score WHERE tag != 'RAW' AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM miranda_score WHERE tag = 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY score desc;");
		}elseif($type == 'inc'){
			$Model = new \Think\Model();
			$ret = $Model->query("SELECT rec_id, gene_symbol FROM miranda_score WHERE tag = 'RAW' AND sig = '{$sig}' AND gene_symbol NOT IN (SELECT gene_symbol FROM miranda_score WHERE tag != 'RAW' AND sig = '{$sig}' GROUP BY gene_symbol) GROUP BY gene_symbol ORDER BY score desc;");
		}
		$this->assign("rlist", $ret);
		$this->assign("sig", $sig);
		$this->assign("miRNA", $ret[0]["mirna"]);
		$this->display("pairList");
	}
	public function geneListS($sig, $type, $p=1){
		if (session("search") == 1){
			$db = M('sms');
		}else{
			$db = M('mirmap_score');
		}
		if ($type == 'raw') {
			$map = array('sig'=>$sig, 'tag'=>'RAW');
		}elseif($type == 'mut'){
			$map = array('sig'=>$sig, 'tag'=>'NEW');
		}
		elseif($type == 'inc'){
			$map = array('sig'=>$sig, 'tag'=>'DIE');
		}
		
		$ret = $db->where($map)->order('gene_symbol asc')->page($_GET['p'].',20')->select();
		$count      = $db->where($map)->count();// 鏌ヨ婊¤冻瑕佹眰鐨勬�昏褰曟暟
		$Page       = new \Think\Page($count,20);// 瀹炰緥鍖栧垎椤电被 浼犲叆鎬昏褰曟暟鍜屾瘡椤垫樉绀虹殑璁板綍鏁�
		$show       = $Page->show();
		$this->assign('page', $show);
		$this->assign("rlist", $ret);
		$this->assign("sig", $sig);
		$this->assign("miRNA", $ret[0]["mirna"]);
		$this->display("pairListS");
	}
}