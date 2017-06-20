object Form5: TForm5
  Left = 282
  Top = 125
  BorderStyle = bsSingle
  Caption = #1055#1086#1089#1090#1088#1086#1077#1085#1080#1077' '#1075#1088#1072#1092#1080#1082#1072
  ClientHeight = 289
  ClientWidth = 555
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -10
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  OnShow = FormShow
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox6: TGroupBox
    Left = 0
    Top = 144
    Width = 209
    Height = 145
    Caption = #1054#1089#1100' Y'
    TabOrder = 0
    object Panel3: TPanel
      Left = 8
      Top = 16
      Width = 185
      Height = 57
      TabOrder = 0
      object Label1: TLabel
        Left = 8
        Top = 8
        Width = 65
        Height = 13
        Caption = #1042#1077#1090#1074#1100' ('#1091#1079#1077#1083'):'
      end
      object comboy1: TComboBox
        Left = 8
        Top = 28
        Width = 169
        Height = 21
        ItemHeight = 0
        TabOrder = 0
        Text = 'comboy1'
        OnChange = comboy1Change
      end
    end
    object Panel4: TPanel
      Left = 8
      Top = 80
      Width = 185
      Height = 57
      TabOrder = 1
      object Label2: TLabel
        Left = 8
        Top = 8
        Width = 54
        Height = 13
        Caption = #1055#1072#1088#1072#1084#1077#1090#1088':'
      end
      object comboy2: TComboBox
        Left = 8
        Top = 28
        Width = 169
        Height = 21
        ItemHeight = 0
        TabOrder = 0
        Text = 'ComboBox2'
      end
    end
  end
  object ListGraph: TListBox
    Left = 216
    Top = 8
    Width = 313
    Height = 193
    ItemHeight = 13
    MultiSelect = True
    TabOrder = 1
  end
  object Button1: TButton
    Left = 216
    Top = 208
    Width = 75
    Height = 25
    Caption = #1044#1086#1073#1072#1074#1080#1090#1100
    TabOrder = 2
    OnClick = Button1Click
  end
  object Button2: TButton
    Left = 456
    Top = 208
    Width = 75
    Height = 25
    Caption = #1043#1088#1072#1092#1080#1082
    TabOrder = 3
    OnClick = Button2Click
  end
  object Button3: TButton
    Left = 376
    Top = 240
    Width = 75
    Height = 25
    Caption = #1054#1090#1084#1077#1085#1072
    TabOrder = 4
    OnClick = Button3Click
  end
  object GroupBox5: TGroupBox
    Left = 0
    Top = 0
    Width = 209
    Height = 145
    Caption = #1054#1089#1100' '#1061
    TabOrder = 5
    object Panel5: TPanel
      Left = 8
      Top = 16
      Width = 185
      Height = 57
      TabOrder = 0
      object Label3: TLabel
        Left = 8
        Top = 8
        Width = 65
        Height = 13
        Caption = #1042#1077#1090#1074#1100' ('#1091#1079#1077#1083'):'
      end
      object combox1: TComboBox
        Left = 8
        Top = 28
        Width = 169
        Height = 21
        ItemHeight = 0
        TabOrder = 0
        Text = 'ComboBox2'
        OnChange = combox1Change
      end
    end
    object Panel6: TPanel
      Left = 8
      Top = 80
      Width = 185
      Height = 57
      TabOrder = 1
      object Label4: TLabel
        Left = 8
        Top = 8
        Width = 54
        Height = 13
        Caption = #1055#1072#1088#1072#1084#1077#1090#1088':'
      end
      object combox2: TComboBox
        Left = 8
        Top = 28
        Width = 169
        Height = 21
        ItemHeight = 0
        TabOrder = 0
        Text = 'ComboBox2'
      end
    end
  end
  object Button4: TButton
    Left = 376
    Top = 208
    Width = 75
    Height = 25
    Caption = #1054#1095#1080#1089#1090#1080#1090#1100' '#1074#1089#1077
    TabOrder = 6
    OnClick = Button4Click
  end
  object Button5: TButton
    Left = 216
    Top = 240
    Width = 75
    Height = 25
    Caption = #1059#1076#1072#1083#1080#1090#1100
    TabOrder = 7
    OnClick = Button5Click
  end
  object Button6: TButton
    Left = 296
    Top = 240
    Width = 75
    Height = 25
    Caption = #1042#1089#1077' '#1089#1074#1103#1079#1080
    TabOrder = 8
    OnClick = Button6Click
  end
  object Button7: TButton
    Left = 296
    Top = 208
    Width = 75
    Height = 25
    Caption = #1042#1089#1077' '#1091#1079#1083#1099
    TabOrder = 9
    OnClick = Button7Click
  end
end
