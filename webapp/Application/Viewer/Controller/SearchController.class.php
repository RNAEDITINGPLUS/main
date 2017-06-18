<?php
namespace Viewer\Controller;
use Think\Controller;
class SearchController extends Controller {
	private $genomicInfo = array();

    public function index(){
    	session('jobId', C('SJOB'));
        $res = $this->engine(I());
		if ($res){
			session('map', $res);
			session('query', I('term'));
			$this->success('ok', U('show'));
		}else{
			$term = I('term');
			$this->error("Your query '{$term}' is not in our database. You can go to GENCODE to validate that {$term} is an official gene symbol or gene. If it is an official symbol and you get this message, then we have not predicted targets for this gene.");
		}
    }
    public function confirm($gs){
    	$query = array('RefSeq_gid'	=>	$gs,);
    	#$res = $this->engine($query);
    	session('query', $gs);
    	session('map', $query);
    	$this->success('Redirecting to result page...', U("show?term=$gs"));
    }
    public function show(){
    	$map = session('map');
    	if ($map){
    		$eNum = 0;
    		$max = 0;
    		$sites = M(C('SCANDIDATES'));
    		//$ret = $sites->where($map)->select();
    		if (isset($map['two_round']) && $map['two_round'] == 1){
				$ret = $sites->where($map['acc'])->select();
				if (count($ret) == 0){
					$ret = $sites->where($map['fuzzy'])->select();
				}
			}else{
				$ret = $sites->where($map)->select();
			}
    		$min = $ret[0]['position'];
    		$chr = str_replace('chr', '', $ret[0]['region']);
    		foreach ($ret as $value) {
    			if ($value['position'] > $max){
    				$max = $value['position'];
    			}elseif ($value['position'] < $min){
    				$min = $value['position'];
    			}
    		}
    		if ($min == $max){
    			$max += 1000;
    			$min -= 1000;
    		}
    		$tNum = count($ret);
    		$result = array(
    			'i'	=>	array(),
    			'u'	=>	array(),
    			'mi'	=>	array(),
    			'mis'	=>	array(),
    		);

    		if (count($ret) < C('MAX_SEARCH_AMOUNT')){
    			$idb = M(C('SINTRON'));
    			$udb = M(C('SUTR'));
    			$midb= M(C('SMIR'));
    			$misdb = M(C('SMIS'));
	    		foreach ($ret as $hit) {
	    			if (strpos($hit['refseq_feat'], 'intron')!==false){
	    				//search intron database
	    				$condition = array(
	    					'job'		=>	$hit['job'],
	    					'pos'		=>	$hit['position'],
	    					'gene'		=>	$hit['refseq_gid'],
	    				);
	    				$t = $idb->where($condition)->select();
	    				$eNum += count($t);
	    				$result['i'] = array_merge($result['i'], $t);
	    			}
	    			if (strpos($hit['refseq_feat'], '3UTR') !== false){
	    				$condition = array(
	    						'job'			=>	$hit['job'],
	    						'edit_pos_raw'		=>	$hit['position'],
	    						'gene'	=>	$hit['refseq_gid'],
	    				);
	    				$t = $udb->where($condition)->select();
	    				$eNum += count($t);
	    				$result['u'] = array_merge($result['u'], $t);
	    			}
	    			if (strpos($hit['refseq_gid'], 'MIR') !== false){
	    				$condition = array(
	    						'job'			=>	$hit['job'],
	    						'edit_pos_raw'		=>	$hit['position'],
	    				);
	    				$t = $midb->where($condition)->select();
	    				foreach ($t as $key => $value) {
	    					if ($value['old_t'] != $value['new_t'] || $value['old_t'] != $value['com_t']){
	    						$eNum += 1;
	    					}
	    				}
	    				$result['mi'] = array_merge($result['mi'], $t);
	    			}
	    			if (strpos($hit['refseq_feat'], 'CDS') !== false){
	    				$condition = array(
	    						'job'			=>	$hit['job'],
	    						'position'		=>	$hit['position'],
	    						'gene'	=>	$hit['refseq_gid'],
	    				);
	    				$t = $misdb->where($condition)->select();
	    				$eNum += count($t);
	    				$result['mis'] = array_merge($result['mis'], $t);
	    			}
	    		}
	    		$this->assign('sites', $ret);
	    		$this->assign('num', $tNum);
	    		$this->assign('enum', $eNum);
	    		$this->assign('up', $max);
	    		$this->assign('dn', $min);
	    		$this->assign('ch', $chr);
	    		$this->assign('ss', $result['i']);
	    		$this->assign('ml', $result['mi']);
	    		$this->assign('ul', $result['u']);
	    		$this->assign('mis', $result['mis']);
	    		$this->display();
    		}else{
    			$this->error('Your query contains more than '.C('MAX_SEARCH_AMOUNT').' results, which has been blocked by REP. Please narrow your search scope.');
    		}
    	}else{

    	}

    }
	public function engine($param){
		$terms = explode(" ", $param['term']);
		$map = array();
		//var_dump(preg_match('/chr.*?\:\d+-\d+/', $param['term']));
		//$this->error();
		if (count($terms)==1){ //god inputs a position combination
			if (preg_match('/chr.*?\:\d+$/', $param['term'])) {
				$tmp = explode(':', $param['term']);
				$map = array(
						'Region'	=> $tmp[0],
						'Position'		=> str_replace('chr', '', $tmp[1]),
				);
			}elseif (preg_match('/chr.*?\:\d+-\d+/', $param['term'])){
				$tmp = explode(':', $param['term']);
				$ran = explode('-', $tmp[1]);
				if (intval($ran[1])-intval($ran[0]) < 0 || intval($ran[1])-intval($ran[0]) > 1000){
					$this->error('Search region should not be greater than 1,000bp.');
				}else{
					$map = array(
							'Region'	=>	$tmp[0],
							'Position'	=>	array('between', array(intval($ran[0]), intval($ran[1]))),
					);
				}
				//var_dump($map);
			}elseif (is_numeric($param['term'])){ //god inputs a gene id

				$symbol = $this->getGeneSymbolbyId($param['term']);

				foreach ($symbol as $geneSymbol) {
					$map['RefSeq_gid'] = $geneSymbol;
					$res = $this->seek($map);
				}

			}else{ //god inputs a gene symbol
				$gsm = F('GSM');
				$term = rtrim(ltrim($param['term']));
				if (strpos($term, 'hsa-')!==FALSE){
					$term = strtolower($term);
				}
				if (strpos($term, '-3p')!==FALSE){
					$term = str_replace('-3p', '', $term);
				}
				if (strpos($term, '-5p')!==FALSE){
					$term = str_replace('-5p', '', $term);
				}
				$term = strtolower($term);
				if (isset($gsm[$term])){
					if (count($gsm[$term]) > 1){ //More than one hit
						foreach ($gsm[$term] as $gs) {
							$this->assign('gs', $gsm[$term]);
							$info = $this->fetch('multihits');
							$data = array(
									'status'	=>	2,
									'info'		=>	$info,
									'url'		=>	'',
							);
							$this->ajaxReturn($data);
						}
					}else{
						$map['two_round'] = 1;
						$map['acc'] = array('RefSeq_gid'	=>	$gsm[$term][0]);
						$map['fuzzy'] = array('RefSeq_gid'	=>	array('like', '%'.$gsm[$term][0].'%'));
						# var_dump($gsm[$term][0]);
						# $map = array('RefSeq_gid'	=>	$gsm[$term][0]);
					}
				}else{
					$this->error('No result.');
				}
				//$genomicInfo = $this->getGeneRegion($map['gene']);
				//$this->assign('gi', $genomicInfo);
			}
		}else{
			$this->error('Sorry, but we are not ready for multiple keywords.');
		}

		return $map;
	}

	private function seek($map){
		$rep = M('rep');
		if (isset($map['two_round']) && $map['two_round'] == 1){
			$ret = $rep->where($map['acc'])->select();
			if (count($ret) == 0){
				$ret = $rep->where($map['fuzzy'])->select();
			}
		}else{
			$ret = $rep->where($map)->select();
		}
		return $ret;
	}

	private function enhenceGeneSearch($term){
		$geneInfo = M('gene_info');
		$map['synonyms'] = array('LIKE', array($term.'%', '%|'.$term.'|%', '%|'.$term), 'OR');
		$data = $geneInfo->where($map)->field('symbol, synonyms')->select();
		return $data;
	}

	private function getGeneSymbolbyId($id){
		$geneInfo = M('gene_info');
		$map['gene_id'] = array('EQ', $id);
		$data = $geneInfo->where($map)->field('symbol')->select();
		return $data;
	}

	private function getGeneRegion($gene){
		$geneInfo = M('gene_info');
		$map['symbol'] = array('EQ', $gene);
		$data = $geneInfo->where($map)->field('chromosome, start, end')->find();
		return $data;
	}
	public function reform(){
		$map = array();
		$file = file("symbol2synonyms_map.txt");
		foreach($file as $line){
			$line = str_replace(PHP_EOL, "", $line);
			$items = explode('{rep}', $line);
			$key = strtolower($items[0]);
			$map[$key] = array($items[0]);
			if ($items[1]!='-'){
				$syns = explode('|', $items[1]);
				foreach ($syns as $syn) {
					$key = strtolower($syn);
					if (!isset($map[$key])){
						$map[$key] = array($items[0]);
					}else{
						$map[$key][] = $items[0];
					}
				}
			}else{
				//echo $items[0];
			}

		}
		F('GSM', $map);
		//print_r($map);
	}
}
?>
