object FormContactR: TFormContactR
  Left = 192
  Top = 107
  Width = 696
  Height = 480
  Caption = 'Контактное сопротивление'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 0
    Top = 0
    Width = 305
    Height = 393
    Caption = 'Контакт-пластина 1'
    TabOrder = 0
    object GroupBox2: TGroupBox
      Left = 0
      Top = 24
      Width = 297
      Height = 193
      Caption = 'Удельные характеристики'
      TabOrder = 0
      object Label8: TLabel
        Left = 8
        Top = 160
        Width = 159
        Height = 13
        Caption = 'Температура плавления Тпл, К'
      end
      object Label7: TLabel
        Left = 8
        Top = 136
        Width = 137
        Height = 13
        Caption = 'Теплопроводность, Вт/мК:'
      end
      object Label6: TLabel
        Left = 8
        Top = 120
        Width = 104
        Height = 13
        Caption = 'сопротивления, 1/С:'
      end
      object Label5: TLabel
        Left = 8
        Top = 104
        Width = 162
        Height = 13
        Caption = 'Темп. коэфф-т электрического '
      end
      object Label3: TLabel
        Left = 8
        Top = 80
        Width = 126
        Height = 13
        Caption = 'Микротвердость Н, Мпа:'
      end
      object Label2: TLabel
        Left = 8
        Top = 56
        Width = 196
        Height = 13
        Caption = 'Удельное сопротивление, р, 10  Ом м:'
      end
      object Label4: TLabel
        Left = 168
        Top = 48
        Width = 7
        Height = 10
        Caption = '-9'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -8
        Font.Name = 'Times New Roman'
        Font.Style = []
        ParentFont = False
      end
      object Label1: TLabel
        Left = 8
        Top = 24
        Width = 53
        Height = 13
        Caption = 'Материал:'
      end
      object Edit5: TEdit
        Left = 216
        Top = 160
        Width = 65
        Height = 21
        TabOrder = 0
        Text = '0'
      end
      object Edit4: TEdit
        Left = 216
        Top = 136
        Width = 65
        Height = 21
        TabOrder = 1
        Text = '0'
      end
      object Edit3: TEdit
        Left = 216
        Top = 112
        Width = 65
        Height = 21
        TabOrder = 2
        Text = '0'
      end
      object Edit2: TEdit
        Left = 216
        Top = 80
        Width = 65
        Height = 21
        TabOrder = 3
        Text = '0'
      end
      object ComboBox1: TComboBox
        Left = 136
        Top = 24
        Width = 145
        Height = 21
        ItemHeight = 13
        TabOrder = 4
        Items.Strings = (
          'Медь'
          'Алюминий')
      end
      object Edit1: TEdit
        Left = 216
        Top = 56
        Width = 65
        Height = 21
        TabOrder = 5
        Text = '0'
      end
    end
    object GroupBox3: TGroupBox
      Left = 0
      Top = 216
      Width = 297
      Height = 169
      Caption = 'Параметры'
      TabOrder = 1
      object Label9: TLabel
        Left = 8
        Top = 16
        Width = 184
        Height = 13
        Caption = 'Средний шаг пов-ти профиля, h, мм:'
      end
      object Label10: TLabel
        Left = 8
        Top = 40
        Width = 205
        Height = 13
        Caption = 'Номинальная площадь контакта, S, мм:'
      end
      object Label11: TLabel
        Left = 212
        Top = 32
        Width = 4
        Height = 10
        Caption = '2'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -8
        Font.Name = 'Times New Roman'
        Font.Style = []
        ParentFont = False
      end
      object Label12: TLabel
        Left = 8
        Top = 64
        Width = 143
        Height = 13
        Caption = 'Сечение проводника, S, мм:'
      end
      object Label13: TLabel
        Left = 148
        Top = 56
        Width = 4
        Height = 10
        Caption = '2'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -8
        Font.Name = 'Times New Roman'
        Font.Style = []
        ParentFont = False
      end
      object Label14: TLabel
        Left = 8
        Top = 88
        Width = 152
        Height = 13
        Caption = 'Толщина контакт-детали, мм:'
      end
      object Edit6: TEdit
        Left = 216
        Top = 16
        Width = 65
        Height = 21
        TabOrder = 0
        Text = '0'
      end
      object Edit7: TEdit
        Left = 216
        Top = 40
        Width = 65
        Height = 21
        TabOrder = 1
        Text = '0'
      end
      object Edit8: TEdit
        Left = 216
        Top = 64
        Width = 65
        Height = 21
        TabOrder = 2
        Text = '0'
      end
      object Edit9: TEdit
        Left = 216
        Top = 88
        Width = 65
        Height = 21
        TabOrder = 3
        Text = '0'
      end
    end
  end
end
