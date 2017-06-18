<?php
namespace Api\Controller;
use Think\Controller;
class RepController extends Controller {
    public function query($positionCombine){
        if(preg_match("/\\d+:\\d+-\\d+/", $positionCombine)){
        	$rep = M('searchable');
        	$c = $this->positionParser($positionCombine);
			if (strpos($c[0], 'chr') === false){
				$map['Region'] = 'chr'.$c[0];
			}else{
				$map['Region'] = $c[0];
			}
			#$Model = new \Think\Model(); // 实例化一个model对象 没有对应任何数据表
			#var_dump("SELECT * FROM candidates WHERE `Job`=".intval(session('jobId'))." AND `Position` BETWEEN ".$c[1]." AND ".$c[2].";");
			#var_dump($Model->query("SELECT * FROM candidates WHERE `Job`=".intval(session('jobId'))." AND `Position` BETWEEN ".$c[1]." AND ".$c[2].";"));
        	$map['Position'] = array(array('gt', $c[1]), array('lt', $c[2]));
            if(strpos($_SERVER['HTTP_REFERER'], 'www.rnaeditplus.org/index.php/Viewer/Search/show') === false){
                $map['Job'] = session('jobId');
            }
        	$tmp = $rep->where($map)->select();
			#var_dump($tmp);
			$ret = $this->resultBuilder($tmp);
    
        	$this->ajaxReturn($ret);
        	/**/
        }else{
        	$this->error("Illegal format!");
        }
    }

    public function privateQuery($positionCombine){
        if(preg_match("/\\d+:\\d+-\\d+/", $positionCombine)){
            $rep = M('candidates');
            $c = $this->positionParser($positionCombine);
            if (strpos($c[0], 'chr') === false){
                $map['Region'] = 'chr'.$c[0];
            }else{
                $map['Region'] = $c[0];
            }
            #$Model = new \Think\Model(); // 实例化一个model对象 没有对应任何数据表
            #var_dump("SELECT * FROM candidates WHERE `Job`=".intval(session('jobId'))." AND `Position` BETWEEN ".$c[1]." AND ".$c[2].";");
            #var_dump($Model->query("SELECT * FROM candidates WHERE `Job`=".intval(session('jobId'))." AND `Position` BETWEEN ".$c[1]." AND ".$c[2].";"));
            $map['Position'] = array(array('gt', $c[1]), array('lt', $c[2]));
            if(strpos($_SERVER['HTTP_REFERER'], 'www.rnaeditplus.org/index.php/Viewer/Search/show') === false){
                $map['Job'] = session('jobId');
            }
            $tmp = $rep->where($map)->select();
            #var_dump($tmp);
            $ret = $this->resultBuilder($tmp);
    
            $this->ajaxReturn($ret);
            /**/
        }else{
            $this->error("Illegal format!");
        }
    }
    
    private function positionParser($combine) {
    	$res = array();
    	
    	$tmp = explode(":", $combine);
    	$res[] = $tmp[0];
    	$tmp = explode("-", $tmp[1]);
    	$res[] = $tmp[0];
    	$res[] = $tmp[1];
    	
    	return $res;
    }
    
    private function resultBuilder($rawArray) {
    	$ret = array();
    	foreach ($rawArray as $raw) {
    		$mut = explode('>', $raw['type']);
    		$ret[] = array(
    				'alleles'	=>	array($mut[0], $mut[1]),
    				'feature_type'	=>	'variation',
    				'assembly_name'	=>	'GRCh38',
    				'clinical_significance'	=>	'',
    				'end'		=>	$raw['position'],
    				'seq_region_name'	=>	$raw['chromosome'],
    				'consequence_type'	=>	$raw['annot'],
    				'strand'	=>	($raw['strand']=='+') ? 1 : 2,
    				'id'		=>	$raw['rec_id'],
    				'start'		=>	$raw['position'],
    		);
    	}
    	
    	return $ret;
    }
}