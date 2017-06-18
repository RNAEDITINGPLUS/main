<?php
namespace Viewer\Controller;
use Think\Controller;
class GOController extends Controller {
	public function dispatcher($g){
		if ($g == 1){
			if(session('reg')){
				$this->goAnalysis(session('reg'));
			}else{
				$this->error('Nothing to be analysed!');
			}
		}elseif ($g == 2){
			if(session('rel')){
				$this->goAnalysis(session('rel'));
			}else{
				$this->error('Nothing to be analysed!');
			}
		}else{
			$this->error('Wrong parameter');
		}
	}
	public function goAnalysis($geneList){
		$genes = explode(', ', $geneList);
		$dict = F('GS2TR');
		$transcripts = array();
		foreach ($genes as $gene) {
			if ($dict[$gene]){
				$transcripts[] = $dict[$gene];
			}
		}
		if (count($transcripts) > 3000){
			$this->assign('tid', $transcripts);
			$this->display('overload');
		}else{
			set_time_limit(0);
			$client = new \SoapClient('http://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl');
			$au = new \Org\Net\authenticate('yys@rnaeditplus.org');
			$autRet = $client->authenticate($au);
			if ($autRet->return == 'true'){
				$inputIds = join(',', $transcripts);
				$idType = 'ENSEMBL_TRANSCRIPT_ID';
				$listName = 'make_up';
				$listType = 0;
				try{
					$al = new \Org\Net\addList($inputIds, $idType, $listName, $listType);
					$addListResponse = $client->addList($al);
				}catch(Exception $e){
					print $e->getMessage();
					exit();
				}
				$thd = 0.1;
				$ct = 2;
				$cr = new \Org\Net\getChartReport($thd, $ct);
				$chartReport = $client->getChartReport($cr);
				$this->assign('rl', $chartReport->return);
				$this->display('goanalysis');
			}else{
				exit($autRet->return);
			}
		}
		
	}
	public function buildENSTDict() {
		$map = array();
		$file = file('OFFICIAL_GENE_SYMBOL2ENSEMBL_TRANSCRIPT_ID.txt');
		foreach ($file as $line) {
			$line = str_replace("\n", "", $line);
			$items = explode("\t", $line);
			$map[$items[0]] = $items[1];
		}
		F('GS2TR', $map);
	}
}
?>
