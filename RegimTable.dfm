object FormR: TFormR
  Left = 242
  Top = 128
  Width = 801
  Height = 444
  AutoSize = True
  Caption = '�������������� �����'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnActivate = FormActivate
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object PageControl: TPageControl
    Left = 0
    Top = 0
    Width = 673
    Height = 393
    ActivePage = TabParams
    TabOrder = 0
    OnChange = PageControlChange
    object TabSheetNodes: TTabSheet
      Caption = '����'
      object StringGrid1: TStringGrid
        Left = 0
        Top = 0
        Width = 665
        Height = 365
        Align = alClient
        ColCount = 12
        DefaultRowHeight = 15
        FixedCols = 0
        Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goColSizing, goEditing]
        TabOrder = 0
        RowHeights = (
          15
          15
          15
          15
          15)
      end
    end
    object TabSheetLinks: TTabSheet
      Caption = '�����'
      ImageIndex = 1
      object StringGrid2: TStringGrid
        Left = 0
        Top = 0
        Width = 665
        Height = 365
        Align = alClient
        ColCount = 10
        DefaultRowHeight = 15
        FixedCols = 0
        Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing]
        TabOrder = 0
      end
    end
    object TabAdditional: TTabSheet
      Caption = '����-�������������'
      ImageIndex = 2
      object StringGrid3: TStringGrid
        Left = 0
        Top = 0
        Width = 665
        Height = 365
        Align = alClient
        ColCount = 10
        DefaultRowHeight = 15
        FixedCols = 0
        Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goColSizing, goEditing]
        TabOrder = 0
        RowHeights = (
          15
          15
          15
          15
          15)
      end
    end
    object TabParams: TTabSheet
      Caption = '���������'
      ImageIndex = 3
      object Label2: TLabel
        Left = 8
        Top = 0
        Width = 65
        Height = 13
        Caption = '��� �������:'
      end
      object ComboBox1: TComboBox
        Left = 80
        Top = 0
        Width = 209
        Height = 21
        ItemHeight = 13
        TabOrder = 0
        Text = '�������������� ����� - ����� �������'
        OnChange = ComboBox1Change
        Items.Strings = (
          '�������������� ����� - ����� �������'
          '�������������� ����� - ������� Z'
          '�������������� ����� - ������������������ �����'
          '�������������� ����� ���������� �� ������������'
          '�������������� ����� - � ��������������� � ������')
      end
      object GroupCommon: TGroupBox
        Left = 8
        Top = 24
        Width = 281
        Height = 129
        Caption = '�����'
        TabOrder = 1
        object Label3: TLabel
          Left = 8
          Top = 16
          Width = 144
          Height = 13
          Caption = '����. ���������� ��������:'
        end
        object Label4: TLabel
          Left = 8
          Top = 32
          Width = 107
          Height = 13
          Caption = '�������� �������, %:'
        end
        object EditMaxIter: TEdit
          Left = 216
          Top = 8
          Width = 50
          Height = 21
          TabOrder = 0
          Text = '50'
        end
        object EditPrec: TEdit
          Left = 216
          Top = 32
          Width = 50
          Height = 21
          TabOrder = 1
          Text = '0,1'
        end
        object CheckTxtFile: TCheckBox
          Left = 8
          Top = 52
          Width = 261
          Height = 17
          Caption = '������� txt ���� � ������������ �������'
          TabOrder = 2
        end
        object CheckSensor: TCheckBox
          Left = 8
          Top = 68
          Width = 153
          Height = 17
          Caption = '����� ��������� �����'
          TabOrder = 3
        end
        object CheckStarFiles: TCheckBox
          Left = 8
          Top = 84
          Width = 225
          Height = 17
          Caption = '����������� ����� ��� �����-������'
          TabOrder = 4
        end
        object ButtStarPath: TBitBtn
          Left = 240
          Top = 80
          Width = 33
          Height = 20
          TabOrder = 5
          OnClick = ButtStarPathClick
          Glyph.Data = {
            36030000424D3603000000000000360000002800000010000000100000000100
            1800000000000003000000000000000000000000000000000000CACACACACACA
            CACACACACACACACACACACACAC6C6C6C5C7C6DCDFDFE2DADB766666335E5C2C3F
            3FC3BFBFCDCDCDC8C8C8CACACACACACACACACACACACAC8C8C8C7C8C8D3D4D5D5
            D1D2A09A9A5367635292906EBDCE0E2122747070D4D5D5C9C9C9CACACACACACA
            C9C9C9C7C8C8D3D4D5D5D1D2A29D9C56696646848575BAD9A5E9FF5B97A6231F
            1D3E3F40A0A0A0CECECECACACAC9C9C9CCCDCDD5D1D1A29C9C566866467D8577
            C8DC9EEEFF9CD3FFA1D5FC5399A28E8B894143436464649A9A9AC9C8C8CBCCCC
            C1BEBF57696845818577C7DC9EF0FF9EE1FF95CEFA95CFF99AE3FA518B90A7A2
            A13836377D7C7CA5A5A5C5C8C8D7D1D18FABAA72BBC5A7F5FF9DEEFF96DCFA96
            D5FB96D9FC9BE6FE81D3EC526A69E1DCDB2E3B3BA6A6A6D6D6D6C8CACACEC9C9
            78BDBDAFF0FC95EEFA96E7FB97E6FD97DFFD96E0FC9BE3FE88DEF34C6868E1DB
            DB3B7273959494D6D6D6D1CCCDB7C4C379C0C6ACECFE92F4FC98E3FD97E8FD97
            E3FD96E2FB9AEDFD71CDD8667B7AE5DCDC6B9DA95A706FD3CECEE7D3D38AB7B6
            91DCE1A2F2FE95EDFC97F5FD98ECFC94E2FC93E7FCB3FFFF5C93969C9C9CFCF6
            F597B2B543696DC1B8B7DACECE87B8B998E8E89FFDFD97F4FC94F0FC95F0FFA7
            F6FFA7E8F387D0D586C1C3D4D7D9ECEDEAABBBCC488897928C87C8C6C56FBABC
            B2FFFF8EFBFB96FBFFA7F5FFA8EFF387CFD582C9CBBCE3E3ECF3F7ECF7FAF8F3
            F3AAB9C2529BAC828784A9BDBC7AC4C6A8FFFFA2FFFFA8F3F387D4D583C9CBBB
            E3E3F2FBFEFEFFFFF4FAFDDFDBDCADBDBA96C9D888DFFA506F6D51A4A4B1EAEC
            B2F7F88AD7D882C9CBBBE1E3F2FBFEFEFFFFF5FAFED6D8D9B4C0BE8AC8C887E8
            EE9DEEFF79C5E060737251B4B58ACFD274B6B78BA8A8F9FEFFFCFFFFF4FAFDD7
            D9D9B1BFBEA6CBCB96DBD9ACF9FCB1F1FFA4E4E836636EA3A19E94BBBB9DBEBE
            CBC6C6B7AFAFF7FFFFDBDCDDB5C1C0A2CACA99DADA7EC9C98BB9B97EBABA77BB
            BA79AEACA49D9CCFD0D0D6CDCDD4CDCDC8C9C9A0A3A3C3BBBB829F9F71C6C691
            D1D18ABABAB5C1C1D5CDCDC9C9C9C3C6C6C9CCCCD4D4D5C8C8C8}
        end
        object CheckCurrReg: TCheckBox
          Left = 8
          Top = 100
          Width = 241
          Height = 17
          Caption = '���������� ������ �� �������� ������'
          TabOrder = 6
        end
      end
      object GroupPredel: TGroupBox
        Left = 296
        Top = 0
        Width = 361
        Height = 361
        Caption = '���������� �����'
        TabOrder = 2
        object RadioGroupPredelType: TRadioGroup
          Left = 8
          Top = 16
          Width = 353
          Height = 65
          Caption = '����� �������'
          ItemIndex = 0
          Items.Strings = (
            '����������� ������������� ��'
            '�������� ���������� ���������� (�������� �)'
            '���������� ���������� �� �������� ����������')
          TabOrder = 0
        end
        object GroupBox2: TGroupBox
          Left = 8
          Top = 80
          Width = 353
          Height = 121
          Caption = '��������� �������'
          TabOrder = 1
          object Label7: TLabel
            Left = 8
            Top = 100
            Width = 87
            Height = 13
            Caption = '������ ������ S:'
          end
          object Label6: TLabel
            Left = 8
            Top = 36
            Width = 176
            Height = 13
            Caption = '��������� �������� ��������� �:'
          end
          object Label5: TLabel
            Left = 8
            Top = 12
            Width = 284
            Height = 13
            Caption = '��������� �������� ������������ ������������� ��:'
          end
          object Label16: TLabel
            Left = 8
            Top = 72
            Width = 256
            Height = 13
            Caption = '(�� ��������� ����������� ������ Vr ���������)'
          end
          object Button4: TButton
            Left = 104
            Top = 96
            Width = 25
            Height = 17
            Caption = '...'
            TabOrder = 0
            OnClick = Button4Click
          end
          object CheckAutoVr0: TCheckBox
            Left = 8
            Top = 56
            Width = 337
            Height = 17
            Caption = '���� ������ ��������� �������� S (����������� ����������)'
            TabOrder = 1
          end
          object EditKn0: TEdit
            Left = 296
            Top = 12
            Width = 50
            Height = 21
            TabOrder = 2
            Text = '0,8'
          end
          object EditT0: TEdit
            Left = 296
            Top = 36
            Width = 50
            Height = 21
            TabOrder = 3
            Text = '15'
          end
          object CheckBox2: TCheckBox
            Left = 8
            Top = 72
            Width = 313
            Height = 17
            Caption = '���� ������ ��������� �������� S (���� �� �����������)'
            TabOrder = 4
          end
        end
        object GroupBox3: TGroupBox
          Left = 8
          Top = 200
          Width = 353
          Height = 89
          Caption = '�������������'
          TabOrder = 2
          object Label8: TLabel
            Left = 240
            Top = 40
            Width = 38
            Height = 13
            Caption = '��� dT:'
          end
          object Label9: TLabel
            Left = 240
            Top = 16
            Width = 44
            Height = 13
            Caption = '��� dKn:'
          end
          object Label10: TLabel
            Left = 176
            Top = 64
            Width = 163
            Height = 13
            Caption = '�� ��������� dT(dKn)=T(Kn)/50'
          end
          object CheckProm: TCheckBox
            Left = 8
            Top = 16
            Width = 209
            Height = 17
            Caption = '���������� ������������� ������'
            TabOrder = 0
            OnClick = CheckPromClick
          end
          object EditdT: TEdit
            Left = 295
            Top = 40
            Width = 50
            Height = 21
            TabOrder = 1
          end
          object CheckGraph: TCheckBox
            Left = 8
            Top = 48
            Width = 145
            Height = 17
            Caption = '������� ���� �������'
            Enabled = False
            TabOrder = 2
          end
          object CheckGrad: TCheckBox
            Left = 8
            Top = 32
            Width = 161
            Height = 17
            Caption = '������ ���������� ������'
            Enabled = False
            TabOrder = 3
          end
          object EditdKn: TEdit
            Left = 295
            Top = 16
            Width = 50
            Height = 21
            TabOrder = 4
          end
        end
        object CheckParamsToFile: TCheckBox
          Left = 16
          Top = 296
          Width = 273
          Height = 17
          Caption = '�������� �������� ��������� ������� � ����'
          TabOrder = 3
        end
        object ButtPrmsFilePath: TBitBtn
          Left = 288
          Top = 296
          Width = 25
          Height = 20
          TabOrder = 4
          OnClick = ButtPrmsFilePathClick
          Glyph.Data = {
            36030000424D3603000000000000360000002800000010000000100000000100
            1800000000000003000000000000000000000000000000000000CACACACACACA
            CACACACACACACACACACACACAC6C6C6C5C7C6DCDFDFE2DADB766666335E5C2C3F
            3FC3BFBFCDCDCDC8C8C8CACACACACACACACACACACACAC8C8C8C7C8C8D3D4D5D5
            D1D2A09A9A5367635292906EBDCE0E2122747070D4D5D5C9C9C9CACACACACACA
            C9C9C9C7C8C8D3D4D5D5D1D2A29D9C56696646848575BAD9A5E9FF5B97A6231F
            1D3E3F40A0A0A0CECECECACACAC9C9C9CCCDCDD5D1D1A29C9C566866467D8577
            C8DC9EEEFF9CD3FFA1D5FC5399A28E8B894143436464649A9A9AC9C8C8CBCCCC
            C1BEBF57696845818577C7DC9EF0FF9EE1FF95CEFA95CFF99AE3FA518B90A7A2
            A13836377D7C7CA5A5A5C5C8C8D7D1D18FABAA72BBC5A7F5FF9DEEFF96DCFA96
            D5FB96D9FC9BE6FE81D3EC526A69E1DCDB2E3B3BA6A6A6D6D6D6C8CACACEC9C9
            78BDBDAFF0FC95EEFA96E7FB97E6FD97DFFD96E0FC9BE3FE88DEF34C6868E1DB
            DB3B7273959494D6D6D6D1CCCDB7C4C379C0C6ACECFE92F4FC98E3FD97E8FD97
            E3FD96E2FB9AEDFD71CDD8667B7AE5DCDC6B9DA95A706FD3CECEE7D3D38AB7B6
            91DCE1A2F2FE95EDFC97F5FD98ECFC94E2FC93E7FCB3FFFF5C93969C9C9CFCF6
            F597B2B543696DC1B8B7DACECE87B8B998E8E89FFDFD97F4FC94F0FC95F0FFA7
            F6FFA7E8F387D0D586C1C3D4D7D9ECEDEAABBBCC488897928C87C8C6C56FBABC
            B2FFFF8EFBFB96FBFFA7F5FFA8EFF387CFD582C9CBBCE3E3ECF3F7ECF7FAF8F3
            F3AAB9C2529BAC828784A9BDBC7AC4C6A8FFFFA2FFFFA8F3F387D4D583C9CBBB
            E3E3F2FBFEFEFFFFF4FAFDDFDBDCADBDBA96C9D888DFFA506F6D51A4A4B1EAEC
            B2F7F88AD7D882C9CBBBE1E3F2FBFEFEFFFFF5FAFED6D8D9B4C0BE8AC8C887E8
            EE9DEEFF79C5E060737251B4B58ACFD274B6B78BA8A8F9FEFFFCFFFFF4FAFDD7
            D9D9B1BFBEA6CBCB96DBD9ACF9FCB1F1FFA4E4E836636EA3A19E94BBBB9DBEBE
            CBC6C6B7AFAFF7FFFFDBDCDDB5C1C0A2CACA99DADA7EC9C98BB9B97EBABA77BB
            BA79AEACA49D9CCFD0D0D6CDCDD4CDCDC8C9C9A0A3A3C3BBBB829F9F71C6C691
            D1D18ABABAB5C1C1D5CDCDC9C9C9C3C6C6C9CCCCD4D4D5C8C8C8}
        end
      end
      object GroupLimits: TGroupBox
        Left = 8
        Top = 152
        Width = 281
        Height = 209
        Caption = '�����������'
        TabOrder = 3
        object CheckQminQmax: TCheckBox
          Left = 8
          Top = 16
          Width = 233
          Height = 17
          Caption = '��������� ����������� Qmin � Qmax'
          TabOrder = 0
        end
        object GroupNebal: TGroupBox
          Left = 0
          Top = 72
          Width = 281
          Height = 137
          Caption = '����������� ���������� ����������'
          TabOrder = 1
          object Label11: TLabel
            Left = 8
            Top = 16
            Width = 168
            Height = 13
            Caption = '�� ����������, dUmax, % U���:'
          end
          object Label12: TLabel
            Left = 8
            Top = 40
            Width = 109
            Height = 13
            Caption = '�� ����, dfmax, ����:'
          end
          object Label13: TLabel
            Left = 8
            Top = 64
            Width = 196
            Height = 13
            Caption = '�� ����-��. �������������,  d��, o.e.:'
          end
          object Label14: TLabel
            Left = 8
            Top = 84
            Width = 98
            Height = 13
            Caption = '�� ��������� �, %:'
          end
          object Label1: TLabel
            Left = 8
            Top = 108
            Width = 177
            Height = 13
            Caption = '�� ���������� ������. �������, %'
          end
          object EditdUmax: TEdit
            Left = 216
            Top = 16
            Width = 50
            Height = 21
            TabOrder = 0
            Text = '5'
          end
          object Editdfmax: TEdit
            Left = 216
            Top = 40
            Width = 50
            Height = 21
            TabOrder = 1
            Text = '5'
          end
          object EditdKnmax: TEdit
            Left = 216
            Top = 60
            Width = 50
            Height = 21
            TabOrder = 2
            Text = '0,1'
          end
          object EditdTmax: TEdit
            Left = 216
            Top = 84
            Width = 50
            Height = 21
            TabOrder = 3
            Text = '30'
          end
          object EditdVrimax: TEdit
            Left = 216
            Top = 108
            Width = 50
            Height = 21
            TabOrder = 4
            Text = '30'
          end
        end
        object CheckNebal: TCheckBox
          Left = 8
          Top = 48
          Width = 233
          Height = 17
          Caption = '������������ ���������:'
          TabOrder = 2
          OnClick = CheckNebalClick
        end
        object CheckBkUse: TCheckBox
          Left = 8
          Top = 32
          Width = 265
          Height = 17
          Caption = '����-������ ����������� ���� (�����-�������)'
          TabOrder = 3
          OnClick = CheckNebalClick
        end
      end
    end
    object TabRusults: TTabSheet
      Caption = '����������'
      ImageIndex = 4
      object Label17: TLabel
        Left = 248
        Top = 24
        Width = 155
        Height = 13
        Caption = '����. �������� ������, �����:'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clGreen
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
      end
      object GroupResCommon: TGroupBox
        Left = 0
        Top = 0
        Width = 201
        Height = 65
        Caption = '�����'
        TabOrder = 0
        object LabelRegim: TLabel
          Left = 8
          Top = 17
          Width = 157
          Height = 16
          Caption = '����� �� ���������!'
          Color = clBtnFace
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clRed
          Font.Height = -13
          Font.Name = 'MS Sans Serif'
          Font.Style = [fsBold]
          ParentColor = False
          ParentFont = False
        end
        object LabelIters: TLabel
          Left = 123
          Top = 40
          Width = 6
          Height = 13
          Caption = '0'
        end
        object LabelItersCaption: TLabel
          Left = 8
          Top = 40
          Width = 112
          Height = 13
          Caption = '���������� ��������:'
        end
      end
      object GroupBox1: TGroupBox
        Left = 0
        Top = 64
        Width = 201
        Height = 153
        Caption = '��������� ����'
        TabOrder = 1
        object Label15: TLabel
          Left = 8
          Top = 16
          Width = 125
          Height = 13
          Caption = '����� ���������� ����:'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object LabelResSensor: TLabel
          Left = 139
          Top = 16
          Width = 6
          Height = 13
          Caption = '0'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object Label18: TLabel
          Left = 8
          Top = 32
          Width = 168
          Height = 13
          Caption = '������, �������������� �����.'
        end
        object Label19: TLabel
          Left = 8
          Top = 48
          Width = 126
          Height = 13
          Caption = '������������ ��������:'
        end
        object StringGridSensor: TStringGrid
          Left = 0
          Top = 64
          Width = 193
          Height = 81
          ColCount = 2
          DefaultRowHeight = 14
          FixedCols = 0
          TabOrder = 0
          ColWidths = (
            81
            87)
        end
      end
      object GroupBox4: TGroupBox
        Left = 208
        Top = 0
        Width = 457
        Height = 217
        Caption = '���������� �����'
        TabOrder = 2
        object LabelResKn: TLabel
          Left = 187
          Top = 24
          Width = 6
          Height = 13
          Caption = '0'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object LabelTpred: TLabel
          Left = 8
          Top = 40
          Width = 88
          Height = 13
          Caption = '�������� �����:'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object LabelResTpred: TLabel
          Left = 187
          Top = 40
          Width = 6
          Height = 13
          Caption = '0'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object Label24: TLabel
          Left = 32
          Top = 56
          Width = 174
          Height = 13
          Caption = '�������� ����������� ��������:'
        end
        object Label25: TLabel
          Left = 272
          Top = 56
          Width = 151
          Height = 13
          Caption = '��������� ������ �� ������:'
        end
        object Label26: TLabel
          Left = 256
          Top = 24
          Width = 155
          Height = 13
          Caption = '����. �������� ������, �����:'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object LabelResGradLink: TLabel
          Left = 417
          Top = 24
          Width = 6
          Height = 13
          Caption = '0'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object LabelKn: TLabel
          Left = 8
          Top = 24
          Width = 171
          Height = 13
          Caption = '����������� ������������� ��:'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clGreen
          Font.Height = -11
          Font.Name = 'MS Sans Serif'
          Font.Style = []
          ParentFont = False
        end
        object StringGridVr: TStringGrid
          Left = 16
          Top = 72
          Width = 233
          Height = 137
          ColCount = 3
          DefaultRowHeight = 14
          FixedCols = 0
          TabOrder = 0
          ColWidths = (
            53
            79
            75)
        end
        object StringGridGradLink: TStringGrid
          Left = 256
          Top = 72
          Width = 193
          Height = 137
          ColCount = 2
          DefaultRowHeight = 14
          FixedCols = 0
          TabOrder = 1
          ColWidths = (
            82
            87)
        end
      end
      object GroupBox5: TGroupBox
        Left = 0
        Top = 216
        Width = 665
        Height = 153
        Caption = '�������� ��������� ���� ������������� ��������:'
        TabOrder = 3
        object StringGrid6: TStringGrid
          Left = 0
          Top = 16
          Width = 665
          Height = 137
          ColCount = 3
          DefaultRowHeight = 14
          FixedCols = 0
          TabOrder = 0
          ColWidths = (
            53
            79
            75)
        end
      end
    end
    object TabRusultsNodes: TTabSheet
      Caption = '����������-����'
      ImageIndex = 5
      object StringGrid5: TStringGrid
        Left = 0
        Top = 0
        Width = 665
        Height = 365
        Align = alClient
        ColCount = 11
        DefaultRowHeight = 15
        Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing]
        TabOrder = 0
      end
    end
    object TabRusultsLinks: TTabSheet
      Caption = '����������-�����'
      ImageIndex = 6
      object StringGrid4: TStringGrid
        Left = 0
        Top = 0
        Width = 665
        Height = 365
        Align = alClient
        ColCount = 20
        DefaultRowHeight = 15
        FixedCols = 4
        Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goColSizing, goEditing]
        TabOrder = 0
        RowHeights = (
          15
          15
          15
          15
          15)
      end
    end
  end
  object PanelButtons: TPanel
    Left = 672
    Top = 0
    Width = 121
    Height = 393
    TabOrder = 1
    object Button2: TButton
      Left = 8
      Top = 24
      Width = 105
      Height = 25
      Caption = '�������� �������'
      TabOrder = 0
      OnClick = Button2Click
    end
    object ButtonRaschet: TButton
      Left = 8
      Top = 120
      Width = 105
      Height = 25
      Caption = '������'
      TabOrder = 1
      OnClick = ButtonRaschetClick
    end
    object Button1: TButton
      Left = 8
      Top = 88
      Width = 105
      Height = 25
      Caption = '��������� �������'
      Enabled = False
      TabOrder = 2
      OnClick = Button1Click
    end
    object Button3: TButton
      Left = 8
      Top = 56
      Width = 105
      Height = 25
      Caption = '������� �������'
      Enabled = False
      TabOrder = 3
      OnClick = Button3Click
    end
    object Button5: TButton
      Left = 8
      Top = 152
      Width = 75
      Height = 25
      Caption = '����� 63'
      TabOrder = 4
      OnClick = Button5Click
    end
    object Button6: TButton
      Left = 8
      Top = 192
      Width = 75
      Height = 25
      Caption = 'Button6'
      TabOrder = 5
      OnClick = Button6Click
    end
    object Button7: TButton
      Left = 8
      Top = 224
      Width = 75
      Height = 25
      Caption = 'Button7'
      TabOrder = 6
      OnClick = Button7Click
    end
  end
  object NewB: TButton
    Left = 8
    Top = 392
    Width = 89
    Height = 25
    Caption = '����� ����'
    Enabled = False
    TabOrder = 2
    OnClick = NewBClick
  end
  object DelB: TButton
    Left = 104
    Top = 392
    Width = 89
    Height = 25
    Caption = '������� ����'
    Enabled = False
    TabOrder = 3
    OnClick = DelBClick
  end
end