# 用户提供SMILES与PubChem CID209254核验

日期：2026-10-24。范围：结构与来源核对，未运行药物靶点预测。

## 结论

用户给出的文本含一个紧接分支括号前的反斜杠：`C[C@@H]\(`。若把该反斜杠按字面写入SMILES，RDKit解析失败；它很可能只是聊天Markdown转义。只移除这一处`\(`中的反斜杠后，得到的SMILES与本轮重新取得的PubChem官方属性API原文**逐字相同**。

```text
C[C@@H]([C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@]5(C4)C6=CC=CC=C6C(=O)N5C)C)C)N(C)C
```

所有@/@@标记保持，未任意赋值或翻转。RDKit解析、化学价检查与立体分配均通过；9个可识别手性中心全部已指定，没有未分配中心；形式电荷0、单一连通分子，分子式C31H46N2O。

生成InChI与官方记录相同，InChIKey均为`SLGWGPQWJRVPAD-PZOAWPBASA-N`。回绘并生成2D MOL/SDF后，再读回MOL的InChIKey仍相同。不同软件的规范SMILES字符顺序可不同，不用字符串形状判断其结构是否改变。

## 能确认与尚不能确认的内容

- 已确认：这是PubChem CID209254所定义的、立体指定完整的结构，可以明确保存/交换。
- 官方当前同义名首项为Spiropachysine，含CAS19587-41-8等；此次取得的完整同义名列表未列精确名称Spiropachysine A。
- 同义名中没有A不证明其与论文药物不同，但该记录本身不能独立完成其与Fang2022实验样品的对应。先前用户要求等待全部结构核准，本轮未将数据库匹配自动登记为研究样品身份已核准。
- 本轮只核验这条字符串与CID；不宣称检验过实物样品或补齐原论文缺失的鉴定证据。后续若用户明确采用CID209254定义研究输入，应记录这是采用的具体结构及对应依据，避免单靠近似名称。

## 文件

- verified_input_smiles.smi：保留用户/官方原字序的正确输入。
- rdkit_canonical_isomeric_smiles.smi：RDKit规范化后的等价异构SMILES。
- PubChem_CID209254_verified_2D.mol及.sdf：含立体信息的二维结构文件。
- CID209254_stereochemistry.png及.svg：含R/S标记的回绘图，已查看。
- PubChem_properties.json、PubChem_synonyms.json：本轮从官方API取得的原始响应。
- Verification_report.json：逐项检验、访问来源、哈希、CIP列表和结论。CIP表的索引是RDKit零起始原子索引，不是天然产物的传统碳编号。

初次写MOL遇到RDKit原生Windows中文路径输出问题，已改为内存生成MolBlock/SDF、由Python写入Unicode路径；未更改结构。原字符串解析失败为本次格式核查的预期结果，不是移除格式转义后的结构失败。
