from pathlib import Path
import sys,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(sys.argv[1]);out=root/'outputs/MFSD10_稀疏与聚合探索_2026-10-26'
old=pd.read_csv(out/'results/Existing_models_reused_all20.tsv',sep='\t');main=pd.read_csv(out/'results/Aggregation_main_all10.tsv',sep='\t');sens=pd.read_csv(out/'results/Aggregation_all60_fixed_specifications.tsv',sep='\t');matched=pd.read_csv(out/'results/Matched_original_singlecell_reused.tsv',sep='\t');b=pd.read_csv(out/'results/Patient_between_parsimonious_sensitivities.tsv',sep='\t');availability=pd.read_csv(out/'results/All_patient_unit_availability.tsv',sep='\t')
rename={'GSE161529':'discovery','GSE176078':'validation'};old['cohort_short']=old.cohort.map(rename)
programs=['MP01','MP02','MP03','MP04','MP05'];cohorts=['discovery','validation']
long=[];wide=[];table=[]
for coh in cohorts:
    for mp in programs:
        sc=old[(old.cohort_short==coh)&(old.program_id==mp)&(old.association=='within')].iloc[0]
        pb=old[(old.cohort_short==coh)&(old.program_id==mp)&(old.association=='between')].iloc[0]
        ag=main[(main.cohort==coh)&(main.program_id==mp)].iloc[0]
        s=sens[(sens.cohort==coh)&(sens.program_id==mp)].set_index('specification')
        r=dict(cohort=coh,program_id=mp)
        for label,z,q in [('original_singlecell',sc,'BH_q'),('patient_pseudobulk',pb,'BH_q'),('within_patient_aggregate',ag,'FDR_main10')]:
            x=dict(cohort=coh,program_id=mp,analysis=label,n_patients=int(z.n_patients),effect_rho=z.effect_rho,ci_low=z.ci_low,ci_high=z.ci_high,p=z.p,FDR=z[q],FDR_family='Old: 5 programs x within/between, per cohort' if label!='within_patient_aggregate' else 'New: 5 programs x 2 cohorts, aggregate main only',status='Reused_without_refitting' if label!='within_patient_aggregate' else 'Exploratory_new_aggregation')
            long.append(x)
            for col in ['n_patients','effect_rho','ci_low','ci_high','p','FDR']:r[label+'_'+col]=x[col]
        raw=s.loc['state_unadjusted','effect_rho'];full=ag.effect_rho
        attenuation=1-abs(full)/abs(raw) if abs(raw)>1e-12 else np.nan
        r['unadjusted_aggregate_rho']=raw;r['depth_only_aggregate_rho']=s.loc['state_depth','effect_rho'];r['depth_mito_without_features_rho']=s.loc['state_depth_mt','effect_rho'];r['original_cellscore_mean_rho']=s.loc['state_cellmean_full','effect_rho'];r['random_same_capacity_rho']=s.loc['random_full','effect_rho'];r['technical_attenuation_fraction']=attenuation
        tags=[]
        if abs(full)<.1:tags.append('效应很弱')
        if ag.ci_low<=0<=ag.ci_high:tags.append('不精确')
        if np.sign(raw)!=np.sign(full):tags.append('技术调整翻向；近零时不作实质解释')
        if abs(raw)>=.1 and attenuation>=.5:tags.append('技术/复杂度调整解释较多（非因果）')
        other=main[(main.cohort!=coh)&(main.program_id==mp)].iloc[0]
        tags.append('聚合跨队列同向' if np.sign(full)==np.sign(other.effect_rho) else '聚合跨队列方向不一致')
        r['interpretation']='；'.join(tags);wide.append(r)
        def fmt(z,q):return f"{z.effect_rho:.3f} [{z.ci_low:.3f},{z.ci_high:.3f}]; q={z[q]:.3g}; n={int(z.n_patients)}"
        label='发现' if coh=='discovery' else '验证'
        technical=f"未调→全调 {raw:.3f}→{full:.3f}；细胞均分 {r['original_cellscore_mean_rho']:.3f}；随机 {r['random_same_capacity_rho']:.3f}"
        table.append(f"| {label} {mp} | {fmt(sc,'BH_q')} | {fmt(pb,'BH_q')} | {fmt(ag,'FDR_main10')} | {technical} |")
pd.DataFrame(long).to_csv(out/'results/Comparison_all30_effects_CI_P_FDR_N.tsv',sep='\t',index=False)
pd.DataFrame(wide).to_csv(out/'results/Comparison_wide_numeric_and_technical_sensitivity.tsv',sep='\t',index=False)
header=['# 三层关联对照表（探索性）','','每格为rho [95%CI]；q；患者数。患者级列为患者间关系，其他两列为患者内关系，不能当作同一估计量。原模型全部复用。旧q为每队列10项BH；聚合主q为两队列10项BH，全部60规格总BH另存。所有原始P见完整30行数值表。','',
'| 队列/程序 | 原单细胞患者内 | 患者级pseudobulk（已做，复用） | 患者内状态聚合 | 技术/评分/随机对照rho |',
'|---|---|---|---|---|',*table,'',
'[完整效应、CI、原始P、FDR、患者数](results/Comparison_all30_effects_CI_P_FDR_N.tsv)｜[全60个聚合规格](results/Aggregation_all60_fixed_specifications.tsv)｜[精简患者间模型](results/Patient_between_parsimonious_sensitivities.tsv)']
(out/'简洁对照表.md').write_text('\n'.join(header),encoding='utf-8')

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
def save(fig,name):
    fig.savefig(out/'figures'/f'{name}.png',dpi=170,bbox_inches='tight');fig.savefig(out/'figures'/f'{name}.svg',bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
for ax,kind,title in zip(axes,['original_singlecell','patient_pseudobulk','within_patient_aggregate'],['Original within-patient single-cell','Between-patient pseudobulk (reused)','Within-patient aggregate (exploratory)']):
    for coh,color,delta,label in [('discovery','#367CA2',-.12,'GSE161529'),('validation','#C17638',.12,'GSE176078')]:
        z=pd.DataFrame(long);z=z[(z.cohort==coh)&(z.analysis==kind)].set_index('program_id').loc[programs]
        for i,(_,r) in enumerate(z.iterrows()):ax.plot([r.ci_low,r.ci_high],[i+delta,i+delta],color=color,lw=1.5)
        ax.scatter(z.effect_rho,np.arange(5)+delta,color=color,s=30,label=label)
    ax.set_yticks(range(5),programs);ax.axvline(0,color='grey',lw=.8);ax.set_title(title,fontsize=10);ax.set_xlabel('Patient-level rho and 95% CI');ax.invert_yaxis();ax.set_xlim((-.1,.13) if kind=='original_singlecell' else (-1,1) if kind=='patient_pseudobulk' else (-.4,.5))
axes[0].legend(fontsize=8,loc='lower left');save(fig,'01_three_level_effect_comparison')
fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained');zero=[]
for ax,coh in zip(axes,cohorts):
    a=availability[(availability.cohort==coh)&availability.primary_patient]
    for _,r in a.iterrows():ax.plot([0,1],[r.singlecell_zero_fraction,r.aggregate_zero_fraction],alpha=.4,color='#367CA2');ax.scatter([0,1],[r.singlecell_zero_fraction,r.aggregate_zero_fraction],s=16,color='#367CA2')
    ax.set_xticks([0,1],['Single cell','Raw-count aggregate']);ax.set_ylim(-.03,1);ax.set_ylabel('Fraction of MFSD10 zeros per patient');ax.set_title(f'{coh}: original primary patients n={len(a)}')
    zero.append(dict(cohort=coh,n_original_patients=len(a),median_singlecell_zero=float(a.singlecell_zero_fraction.median()),median_aggregate_zero=float(a.aggregate_zero_fraction.median()),n_aggregate_effect_patients=int(main[main.cohort==coh].n_patients.iloc[0])))
fig.suptitle('Zero reduction is a measurement diagnostic, not proof of association');save(fig,'02_zero_fraction_reduction');pd.DataFrame(zero).to_csv(out/'results/Zero_fraction_diagnostics.tsv',sep='\t',index=False)
specorder=['state_unadjusted','state_depth','state_depth_mt','state_full','state_cellmean_full','random_full'];labels=['No adjustment','Depth only','Depth + mito','Full adjustment','Mean original score','Random units']
fig,axes=plt.subplots(1,2,figsize=(13,4.5),layout='constrained')
for ax,coh in zip(axes,cohorts):
    z=sens[sens.cohort==coh].pivot(index='program_id',columns='specification',values='effect_rho').loc[programs,specorder]
    im=ax.imshow(z,vmin=-.5,vmax=.5,cmap='RdBu_r',aspect='auto');ax.set_xticks(range(6),labels,rotation=35,ha='right',fontsize=8);ax.set_yticks(range(5),programs);ax.set_title(coh)
    for i,row in enumerate(z.to_numpy()):
        for j,v in enumerate(row):ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=8)
    fig.colorbar(im,ax=ax,shrink=.7)
save(fig,'03_fixed_technical_sensitivities')

allp=pd.read_csv(out/'results/Aggregation_patient_effects_and_failures.tsv',sep='\t');failed=allp[(allp.specification=='state_full')&allp.primary_patient & allp.reason_not_estimated.notna()].drop_duplicates(['cohort','patient_id'])
failed.to_csv(out/'results/Primary_patients_not_estimable.tsv',sep='\t',index=False)
mp02=main[main.program_id=='MP02'].set_index('cohort');mp05=main[main.program_id=='MP05'].set_index('cohort')
stop={'analysis_type':'exploratory_reanalysis','stop_automatic_expression_program_association_optimization':True,'reason':'Zero fractions fall and validation MP02/MP05 strengthen, but discovery estimates remain imprecise and all five discovery aggregate tests fail multiplicity correction; patient-between models remain unsupported. No robust cross-cohort association established.','NMF_rerun':False,'signature_changed':False,'granularity_changed_after_effects':False,'primary_patient_counts':{'discovery':13,'validation_original':16,'validation_aggregate':15},'no_additional_analyses_authorized_by_this_result':True}
(out/'audit/STOP_THIS_ROUTE.json').write_text(json.dumps(stop,indent=2),encoding='utf-8')
report=['# MFSD10稀疏与聚合探索性再分析','',
'**结论：聚合明显减少了观测零值，并增强了验证队列MP02/MP05的患者内关联；但没有补齐发现队列证据，整体尚不足以建立可靠的跨队列MFSD10关联。按本轮停止规则，不再继续优化这条表达程序关联路线。**','',
'两个队列已查看，本轮所有新增结果均为探索性；不是新的预先独立确认性验证。5个程序及旧结果全部保留，未重跑NMF、未筛新signature、未做MAGIC/ALRA补值。','',
'## 既有模型核查与复用','',
'患者间分析原本已做：每患者恶性RNA counts求和、MFSD10 log2(CPM+1)与该患者原固定程序分数均值作偏秩相关。原13/16患者、IHC/深度/线粒体调整、效应CI/P/FDR全部原样复用；本轮全基因pseudobulk仅作原始counts审计，目标标准化数值和患者程序均分逐项与旧表相符。原单细胞患者内模型先控制UMI、基因数、线粒体，再每患者估计效应、患者级汇总，未把细胞当独立患者。','',
'[原20项模型完整列表](results/Existing_models_reused_all20.tsv)｜[固定方案与方法](探索性再分析_固定规则.md)','',
'## 非重叠状态聚合','',
'每患者×样本单独构建约20细胞的容量约束邻域，MFSD10从状态特征、归一化及PCA输入排除，PCA状态另去除深度/线粒体线性成分；只用RNA counts，不用整合矩阵/UMAP。每细胞恰好一个单元，全基因UMI守恒。全部发现singlet形成425单元、验证作者候选形成1215单元；随机同容量分组另存。状态分组的平均PCA平方距离低于随机分组，说明确实聚合了更相近的状态。','',
'冻结50基因的原log1p(CP10K)均分函数应用于单元counts，目标从分母排除；固定注释覆盖94%–100%。逐单元细胞数、总UMI、检出基因数和程序注释/实际检出覆盖均保存。非线性变换意味着聚合counts评分与原细胞均分不是同一个量，因此后者也单列敏感性，不把它们混同。','',
'发现主要13患者均可估计；验证原16Naïve患者中，4465只有124细胞、6单元，按固定不足8单元规则不估计，主要聚合为15患者。其他患者未因MFSD10或P值筛选；不足的原因为单元数或残差自由度，详见完整表。相同13/15患者的旧效应也复用对照。','',
'## 简洁对照','',
'[打开完整三层简洁对照表](简洁对照表.md)｜[30项效应/CI/原始P/FDR/N](results/Comparison_all30_effects_CI_P_FDR_N.tsv)｜[逐患者效应与不可估原因](results/Aggregation_patient_effects_and_failures.tsv)','',
'![三层效应](figures/01_three_level_effect_comparison.png)','',
'| 程序 | 发现聚合rho [95%CI]；q | 验证聚合rho [95%CI]；q | 判断 |','|---|---|---|---|',
*[(lambda p: f"| {p} | "+' | '.join(f"{r.effect_rho:.3f} [{r.ci_low:.3f},{r.ci_high:.3f}]；{r.FDR_main10:.3g}" for _,r in main[main.program_id==p].set_index('cohort').loc[cohorts].iterrows())+' | '+{'MP01':'近零，跨队列方向不一致','MP02':'正向一致，但发现不精确','MP03':'效应很弱，区间跨零','MP04':'发现不精确；对评分/技术规格较敏感','MP05':'正向一致，但发现不精确；外部增强'}[p]+' |')(p) for p in programs],'',
'原单细胞效应大多很弱；患者间点估计可有中等幅度，但CI很宽，不能简单说“效应都小”。精简到IHC+一个深度指标后，两队列患者间关联仍无多重校正支持，没有在13/16患者上追加一长串变量。','',
'## 稀疏与技术因素是否解释了变化','',
f"原主要患者集的患者零值比例中位数：发现{zero[0]['median_singlecell_zero']:.1%}→{zero[0]['median_aggregate_zero']:.1%}，验证{zero[1]['median_singlecell_zero']:.1%}→{zero[1]['median_aggregate_zero']:.1%}。这是观测计数更充分，不等于每个细胞真的表达，也不等于恢复了单细胞真值。",'',
'验证MP02/MP05的rho升至约0.285/0.282，主10项BH均约0.00137；相同患者子集的旧结果比较，增强并非仅由少了一名患者造成。发现则约0.088/0.136、CI跨零、主q约0.379/0.370，不能把验证单队列增强说成跨队列证据已建立。','',
'六个预定规格全部报告：不调整、仅深度、深度+线粒体、完整三技术调整、原细胞分数均值、随机同容量聚合。MP02/MP05在验证中的增强不只依赖一个评分变换或不控制深度；随机聚合信号较弱。不过发现效应仍不精确，无法据此证明稀疏性就是此前不显著的原因。弱信号中的翻向或调整敏感，不能直接断言是技术因果。','',
'没有中等以上未调整效应符合本轮“调整后衰减至少一半”的技术解释较多判据；不为凑类别强行贴此标签。MP01/MP03主要是很弱/不稳定，MP02/MP05是同向但发现不精确，MP04对评分/调整规格敏感。全部60规格及总BH保留，未选择最显著规格作为结论。','',
'![零值](figures/02_zero_fraction_reduction.png)','',
'![固定技术敏感性](figures/03_fixed_technical_sensitivities.png)','',
'[相同患者配对增幅CI](results/Matched_patient_delta_vs_singlecell.tsv)｜[全部60规格](results/Aggregation_all60_fixed_specifications.tsv)｜[患者间精简敏感性](results/Patient_between_parsimonious_sensitivities.tsv)','',
'**停止判断：测量稀疏性改善了，外部队列部分关联也增强了，但没有实质补足两队列共同支持这一关键缺口。到此停止调整聚合粒度、邻域、signature或模型，不自动追加其他分析；这不是取消MFSD10候选方向，也不是证明其无作用。**','',
'复现见scripts/run.ps1；聚合原始counts及患者pseudobulk保存在objects；所有单元成员、状态特征、不可估原因、固定方案及哈希在metadata/audit。旧结果和5个程序均保持。']
(out/'结果摘要_探索性聚合再分析.md').write_text('\n'.join(report),encoding='utf-8')
p=out/'audit/PROGRESS.json';d=json.loads(p.read_text());d.update(comparison_complete=True,stage='final_audit',stop_optimization=True);p.write_text(json.dumps(d,indent=2),encoding='utf-8')
print(json.dumps(stop,indent=2));print('COMPARISON_AND_REPORT_COMPLETE')
